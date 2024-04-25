%% Descend
clc 
clear

%count atom in each direction
N=8;

% Create massive
work_arr = zeros(N, N, N, 6);
% Create corner points's deviation
deviation = 0.3;
work_arr(1, 1, 1, 1) = -deviation;
work_arr(N, 1, 1, 1) = deviation;
work_arr(1, N, 1, 1) = -deviation;
work_arr(N, N, 1, 1) = deviation;
work_arr(1, 1, 1, 2) = -deviation;
work_arr(N, 1, 1, 2) = -deviation;
work_arr(1, N, 1, 2) = deviation;
work_arr(N, N, 1, 2) = deviation;

% Transform 4D massive in 1D (needed? maybe not need, because methods have
% changed)
opt_arr = reshape(work_arr, [1, 6*N*N*N]);

% Calculate distance between atoms
L = 1;
dx = L/N;
dy = L/N;
dz = L/N;

% Constants for elastic-elastic
c11 = 27.5;
c12 = 17.9;
c44 = 5.43;

% Constants for diff polarization-polarization
G11 = 51;
G14 = 0;
G44 = 2;

% Constant for polarization
a1 = -3.712;
a11 = 6.079;
a12 = 1.303;
a111 = 1.294;
a112 = -1.95;
a123 = -2.5;
a1111 = 3.863;
a1112 = 2.529;
a1122 = 1.637;
a1123 = 1.367;

% Constants for elastic-polarization
q11 = 14.2;
q12 = -0.74;
q44 = 1.57;

% functions energy for all system
energy_domen_3d = @(domen_arr) sum( ...
    domen_arr_3d(a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, ...
    a1123, domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), ...
    domen_arr(:, :, :, 3)), "all");

energy_elasticity_3d = @(a_arr) sum(elasticity_arr_3d(c11,c12,c44,...
    func_u_3d(func_arr_3d(a_arr, N), 1, 1, dx, dx), ...
    func_u_3d(func_arr_3d(a_arr, N), 2, 2, dy, dy), ...
    func_u_3d(func_arr_3d(a_arr, N), 3, 3, dz, dz), ...
    func_u_3d(func_arr_3d(a_arr, N), 1, 2, dx, dy), ...
    func_u_3d(func_arr_3d(a_arr, N), 2, 3, dy, dz), ...
    func_u_3d(func_arr_3d(a_arr, N), 1, 3, dx, dz)...
    ), "all");

energy_interaction_3d = @(a_arr, domen_arr) sum(interaction_arr_3d(q11, q12, q44, ...
    func_u_3d(a_arr, 1, 1, dx, dx), ...
    func_u_3d(a_arr, 2, 2, dy, dy), ...
    func_u_3d(a_arr, 3, 3, dz, dz), ...
    func_u_3d(a_arr, 1, 2, dx, dy), ...
    func_u_3d(a_arr, 2, 3, dy, dz), ...
    func_u_3d(a_arr, 1, 3, dx, dz), ...
    domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), domen_arr(:, :, :, 3) ...
    ), "all");

energy_domen_gradient_3d = @(domen_arr) sum(domen_gradient_arr_3d(G11, G14, G44, ...
    domen_diff_3d(domen_arr(:, :, :, 1), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 3, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 1), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 3, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 1), 3, dx) ...
    ), "all");

optfunc = @(a_arr) energy_elasticity_3d(func_arr_3d(a_arr(1:numel(a_arr)/2), N)) ... 
    + energy_domen_3d(func_arr_3d(a_arr(numel(a_arr)/2+1:end), N)) ...
    + energy_interaction_3d(func_arr_3d(a_arr(1:numel(a_arr)/2), N), ...
                            func_arr_3d(a_arr(numel(a_arr)/2+1:end), N)) ...
    + energy_domen_gradient_3d(func_arr_3d(a_arr(numel(a_arr)/2+1:end), N));

% functions energy for system 3х3х3 for gradient
grad_energy_domen_3d = @(domen_arr) sum( ...
    domen_arr_3d(a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, ...
    a1123, domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), ...
    domen_arr(:, :, :, 3)), "all");

grad_energy_elasticity_3d = @(a_arr) sum(elasticity_arr_3d(c11,c12,c44,...
    func_u11_3d(a_arr(:, :, :, 1), dx), ...
    func_u22_3d(a_arr(:, :, :, 2), dy), ...
    func_u33_3d(a_arr(:, :, :, 3), dz), ...
    func_u12_3d(a_arr(:, :, :, 1), a_arr(:, :, :, 2), dx, dy), ...
    func_u23_3d(a_arr(:, :, :, 2), a_arr(:, :, :, 3), dy, dz), ...
    func_u13_3d(a_arr(:, :, :, 1), a_arr(:, :, :, 3), dx, dz)...
    ), "all");

grad_energy_interaction_3d = @(a_arr, domen_arr) sum(interaction_arr_3d(q11, q12, q44, ...
    func_u11_3d(a_arr(:, :, :, 1), dx), ...
    func_u22_3d(a_arr(:, :, :, 2), dy), ...
    func_u33_3d(a_arr(:, :, :, 3), dz), ...
    func_u12_3d(a_arr(:, :, :, 1), a_arr(:, :, :, 2), dx, dy), ...
    func_u23_3d(a_arr(:, :, :, 2), a_arr(:, :, :, 3), dy, dz), ...
    func_u13_3d(a_arr(:, :, :, 1), a_arr(:, :, :, 3), dx, dz), ...
    domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), domen_arr(:, :, :, 3) ...
    ), "all");

grad_energy_domen_gradient_3d = @(domen_arr) sum(domen_gradient_arr_3d(G11, G14, G44, ...
    domen_diff_1_3d(domen_arr(:, :, :, 1), dx), ...
    domen_diff_2_3d(domen_arr(:, :, :, 2), dy), ...
    domen_diff_3_3d(domen_arr(:, :, :, 3), dz), ...
    domen_diff_2_3d(domen_arr(:, :, :, 1), dy), ...
    domen_diff_1_3d(domen_arr(:, :, :, 2), dy), ...
    domen_diff_3_3d(domen_arr(:, :, :, 2), dz), ...
    domen_diff_2_3d(domen_arr(:, :, :, 3), dy), ...
    domen_diff_1_3d(domen_arr(:, :, :, 3), dx), ...
    domen_diff_3_3d(domen_arr(:, :, :, 1), dz) ...
    ), "all");

gradfunc = @(a_arr) grad_energy_elasticity_3d(a_arr(:, :, :, 1:3)) ... 
    + grad_energy_domen_3d(a_arr(:, :, :, 4:6)) ...
    + grad_energy_interaction_3d(a_arr(:, :, :, 1:3), ...
                            a_arr(:, :, :, 4:6)) ...
    + grad_energy_domen_gradient_3d(a_arr(:, :, :, 4:6));

% boundary conditions
pin = 0*work_arr + 1;
% fixetion corner points
for i = 1:3
    pin(1, 1, 1, i) = 0;
    pin(N, 1, 1, i) = 0;
    pin(1, N, 1, i) = 0;
    pin(N, N, 1, i) = 0;
end
% fixetion lower face
for i = 1:3
    for j = 1:3
        pin(i, j, 1, 3) = 0;
    end
end
pin = reshape(pin, [1, N*N*N*6]);

% Enter system's energy and optimization 
optfunc(opt_arr)
my_arr2 = descend_with_pin_and_step(gradfunc, opt_arr, pin, 0.1);
optfunc(my_arr2)
my_arr2 = descend_with_pin_and_step(gradfunc, my_arr2, pin, 0.01);
optfunc(my_arr2)

arr = my_arr2;
%% 
my_arr2 = reshape(my_arr2, [N, N, N, 6]);
%% Draw 3D
draw_arr = my_arr2;

value_arr_x = zeros(N, N, N);
value_arr_y = zeros(N, N, N);
value_arr_z = zeros(N, N, N);
value_arr_domen_x = zeros(N, N, N);
value_arr_domen_y = zeros(N, N, N);
value_arr_domen_z = zeros(N, N, N);
for i = 1:N
    value_arr_x(i, :, :) = draw_arr(i, :, :, 1) + i-1;
    value_arr_y(:, i, :) = draw_arr(:, i, :, 2) + i-1;
    value_arr_z(:, :, i) = draw_arr(:, :, i, 3) + i-1;
    value_arr_domen_x(i, :, :) = draw_arr(i, :, :, 4);
    value_arr_domen_y(:, i, :) = draw_arr(:, i, :, 5);
    value_arr_domen_z(:, :, i) = draw_arr(:, i, :, 6);
end
n = 0;
figure;
hold on;
xlabel('X');
ylabel('Y');
zlabel('Z');

for i = 1:numel(value_arr_x(:, 1, 1))
    for j = 1:numel(value_arr_x(1, :, 1))
        for k = 1:numel(value_arr_x(1, 1, :))
            n = n+1;
            plot3(value_arr_x(i, j, k), value_arr_y(i, j, k), value_arr_z(i, j, k), 'o');
            quiver3(value_arr_x(i, j, k), value_arr_y(i, j, k), value_arr_z(i, j, k), ...
                value_arr_domen_x(i, j, k)*10, value_arr_domen_y(i, j, k)*10, value_arr_domen_z(i, j, k)*10, 0.2);
%             plot3(value_arr_x(i, j, 1), value_arr_y(i, j, 1), value_arr_z(i, j, 1), 'o');
        end
    end
end
hold off;
%% functions
% Transform elastic or domen 1D massive in 4D massive (N, N, N, 3)
function array = func_arr_3d(a_arr, N)
    array = reshape(a_arr, [N, N, N, 3]);
end

function array = func_u_3d(a_arr, index1, index2, d1, d2)
    num_x = numel(a_arr(:, 1, 1, 1));
    num_y = numel(a_arr(1, :, 1, 1));
    num_z = numel(a_arr(1, 1, :, 1));
    array = zeros(num_x, num_y, num_z);
    a = a_arr(:, :, :, index1);
    for i=1:(index2-1)
        a = permute(a, [2 3 1]);
    end
    array(1:(numel(array(:, 1, 1))-1), :, :) = array(1:(numel(array(:, 1, 1))-1), :, :)...
        + diff(a(:, :, :), 1, 1)/d2;
    for i=1:(4-index2)
        array = permute(array, [2 3 1]);
    end
    if index1 ~= index2
        a = a_arr(:, :, :, index2);
        for i=1:(index1-1)
            a = permute(a, [2 3 1]);
            array = permute(array, [2 3 1]);
        end
        array(1:(numel(array(:, 1, 1))-1), :, :) = array(1:(numel(array(:, 1, 1))-1), :, :)...
            + diff(a(:, :, :), 1, 1)/d1;
        for i=1:(4-index1)
            array = permute(array, [2 3 1]);
        end
    end
end

% Calculation strain in different direction (3, 3, 3, 1) massive
% Is that exactly mathematically correct?
% Maybe? Howewer need checking 
% function array = func_u11_3d(arr_x, d_x)
%     array = cat(1, arr_x(2,:,:)-arr_x(1,:,:), arr_x(3,:,:)-arr_x(1,:,:), arr_x(3,:,:)-arr_x(2,:,:))/d_x;
% end
% function array = func_u22_3d(arr_y, d_y)
%     array = cat(2, arr_y(:,2,:)-arr_y(:,1,:), arr_y(:,3,:)-arr_y(:,1,:), arr_y(:,3,:)-arr_y(:,2,:))/d_y;
% end
% function array = func_u33_3d(arr_z, d_z)
%     array = cat(3, arr_z(:,:,2)-arr_z(:,:,1), arr_z(:,:,3)-arr_z(:,:,1), arr_z(:,:,3)-arr_z(:,:,2))/d_z;
% end
% function array = func_u12_3d(arr_x, arr_y, d_x, d_y)
%     array = cat(1, arr_y(2,:,:)-arr_y(1,:,:), arr_y(3,:,:)-arr_y(1,:,:), arr_y(3,:,:)-arr_y(2,:,:))/d_x + ...
%             cat(2, arr_x(:,2,:)-arr_x(:,1,:), arr_x(:,3,:)-arr_x(:,1,:), arr_x(:,3,:)-arr_x(:,2,:))/d_y;
% end
% function array = func_u13_3d(arr_x, arr_z, d_x, d_z)
%      array = cat(1, arr_z(2,:,:)-arr_z(1,:,:), arr_z(3,:,:)-arr_z(1,:,:), arr_z(3,:,:)-arr_z(2,:,:))/d_x + ...
%          cat(3, arr_x(:,:,2)-arr_x(:,:,1), arr_x(:,:,3)-arr_x(:,:,1), arr_x(:,:,3)-arr_x(:,:,2))/d_z;
% end
% function array = func_u23_3d(arr_y, arr_z, d_y, d_z)
%     array = cat(2, arr_z(:,2,:)-arr_z(:,1,:), arr_z(:,3,:)-arr_z(:,1,:), arr_z(:,3,:)-arr_z(:,2,:))/d_y + ...
%          cat(3, arr_y(:,:,2)-arr_y(:,:,1), arr_y(:,:,3)-arr_y(:,:,1), arr_y(:,:,3)-arr_y(:,:,2))/d_z;
% end
% % This could allow vectorization to be applied

% function array = func_u11_3d(arr_x, d_x)
%     array = cat(1, diff(arr_x(:, :, :), 1, 1)/d_x, zeros(1, 3, 3));
% end


function array = func_u11_3d(arr_x, d_x)
    array = zeros(3, 3, 3);
    array(1:2, :, :) = array(1:2, :, :)...
        + diff(arr_x(:, :, :), 1, 1)/d_x;
end


function array = func_u22_3d(arr_y, d_y)
    array = zeros(3, 3, 3);
    array(:, 1:2, :) = array(:, 1:2, :)...
        + diff(arr_y(:, :, :), 1, 2)/d_y;
end


function array = func_u33_3d(arr_z, d_z)
    array = zeros(3, 3, 3);
    array(:, :, 1:2) = array(:, :, 1:2)...
        + diff(arr_z(:, :, :), 1, 3)/d_z;
end


function array = func_u12_3d(arr_x, arr_y, d_x, d_y)
    array = zeros(3, 3, 3);
    array(1:2, :, :) = array(1:2, :, :)...
        + diff(arr_y(:, :, :), 1, 1)/d_x;
    array(:, 1:2, :) = array(:, 1:2, :)...
        + diff(arr_x(:, :, :), 1, 2)/d_y;
end


function array = func_u13_3d(arr_x, arr_z, d_x, d_z)
    array = zeros(3, 3, 3);
    array(1:2, :, :) = array(1:2, :, :)...
        + diff(arr_z(:, :, :), 1, 1)/d_x;
    array(:, :, 1:2) = array(:, :, 1:2)...
        + diff(arr_x(:, :, :), 1, 3)/d_z;
end


function array = func_u23_3d(arr_y, arr_z, d_y, d_z)
    array = zeros(3, 3, 3);
    array(:, 1:2, :) = array(:, 1:2, :)...
        + diff(arr_z(:, :, :), 1, 2)/d_y;
    array(:, :, 1:2) = array(:, :, 1:2)...
        + diff(arr_y(:, :, :), 1, 3)/d_z;
end

% Calculation domen-domen interaction

function array = domen_diff_3d(p, index, d)
    num_x = numel(p(:, 1, 1));
    num_y = numel(p(1, :, 1));
    num_z = numel(p(1, 1, :));
    array = zeros(num_x, num_y, num_z);
    a = p;
    for i=1:(index-1)
        a = permute(a, [2 3 1]);
    end
    array(1:(numel(array(:, 1, 1))-1), :, :) = array(1:(numel(array(:, 1, 1))-1), :, :)...
        + diff(a(:, :, :), 1, 1)/d;
    for i=1:(4-index)
        array = permute(array, [2 3 1]);
    end
end

function array = domen_diff_1_3d(p, d)
    array = zeros(3, 3, 3);
    array(1:2, :, :) = array(1:2, :, :) ...
        + diff(p, 1, 1)/d;
end


function array = domen_diff_2_3d(p, d)
    array = zeros(3, 3, 3);
    array(:, 1:2, :) = array(:, 1:2, :) ...
        + diff(p, 1, 2)/d;
end


function array = domen_diff_3_3d(p, d)
    array = zeros(3, 3, 3);
    array(:, :, 1:2) = array(:, :, 1:2) ...
        + diff(p, 1, 3)/d;
end

% Calculation domen energy
function array = domen_arr_3d(a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, a1123, p1, p2, p3)
    p1_2 = p1.^2;
    p2_2 = p2.^2;
    p3_2 = p3.^2;
    p1_4 = p1_2.^2;
    p2_4 = p2_2.^2;
    p3_4 = p3_2.^2;
    p1_6 = p1_2.^3;
    p2_6 = p2_2.^3;
    p3_6 = p3_2.^3;
    array = a1*(p1_2 + p2_2 + p3_2) ...
          + a11*(p1_4 + p2_4 + p3_4) ...
          + a12*(p1.*p2 + p2.*p3 + p1.*p3) ...
          + a111*(p1_6 + p2_6 + p3_6) ...
          + a112*(p1_4.*(p2_2 + p3_2) + p2_4.*(p1_2 + p3_2) + p3_4.*(p1_2 + p2_2)) ...
          + a123*(p1_2.*p2_2.*p3_2) ...
          + a1111*(p1_4.^2 + p2_4.^2 + p3_4.^2) ...
          + a1112*(p1_6.*(p2_2 + p3_2) + p2_6.*(p1_2 + p3_2) + p3_6.*(p1_2 + p2_2)) ...
          + a1122*(p1_4.*p2_4 + p2_4.*p3_4 + p1_4.*p3_4) ...
          + a1123*(p1_4.*p2_2.*p3_2 + p1_2.*p2_4.*p3_2 + p1_2.*p2_2.*p3_4);
end

% Calculation elastic enegry
function array = elasticity_arr_3d(c11,c12,c44,u11,u22,u33, u12, u23, u13)
    array = 1/2 * c11*(u11.^2 + u22.^2 + u33.^2) ...
            + c12*(u11.*u22 + u22.*u33 + u11.*u33) ...
            + 2*c44*(u12.^2 + u23.^2 + u13.^2);
end

% Calculation elastic-domen energy
function array = interaction_arr_3d(q11, q12, q44, u11, u22, u33, u12, u23, u13, p1, p2, p3)
    array = -q11*(u11.*p1.^2 + u22.*p2.^2 + u33.*p3.^2)...
            -q12*(u11.*(p2.^2+p3.^2) + u22.*(p1.^2+p3.^2) + u33.*(p1.^2+p2.^2))...
            -2*q44*(u12.*p1.*p2 + u23.*p2.*p3 + u13.*p1.*p3);
end

% Calculation domen-domen energy
function array = domen_gradient_arr_3d(G11, G14, G44, p11, p22, p33, p12, p21, p23, p32, p31, p13)
    array = 1/2*G11*(p11.^2 + p22.^2 + p33.^2) ...
            + G14*(p11.*p22 + p22.*p33 + p11.*p33) ...
            + G44*(p12.^2 + p21.^2 + p23.^2 + p32.^2 + p31.^2 + p13.^2);
end

% Gradient descend with boundary conditions (fixetion point)
function x_cur = descend_with_pin_and_step(F, x0, pin, step_mag)
x_cur = x0;
dxmag = 0.00001;
for i_descend=1:20 % Adjust as you like
    g = new_grad(F,x_cur,dxmag);
    g = g.*pin;
    g_dir = g/norm(g);
    dx = -g_dir*step_mag;
    x_cur = x_cur+dx;
    
end
end

% Spatial (?) gradient 
% (need add different step for different parameters)
function g = new_grad(F, x0, dxmag)
n = cast((numel(x0)/6)^(1/3), "int8"); %разбиение
x = reshape(x0, [n, n, n, 6]);
g = zeros(n, n, n, 6);
% consider inside the volume
for i = 2:n-1
    for j = 2:n-1
        for k = 2:n-1
            for l = 1:6
                F0 = F(x(i-1:i+1, j-1:j+1, k-1:k+1, :));
                x_cur = x;
                x_cur(i, j, k, l) = x(i, j, k, l) + dxmag;
                F1 = F(x_cur(i-1:i+1, j-1:j+1, k-1:k+1, :));
                g(i, j, k, l) = (F1-F0)/dxmag;
            end
        end
    end
end
% consider on bound (maybe need transform in one loop)
for j = 2:n-1
    for k = 2:n-1
        for l = 1:6
            F0 = F(x(1:3, j-1:j+1, k-1:k+1, :));
            x_cur = x;
            x_cur(1, j, k, l) = x(1, j, k, l) + dxmag;
            F1 = F(x_cur(1:3, j-1:j+1, k-1:k+1, :));
            g(1, j, k, l) = (F1-F0)/dxmag;

            F0 = F(x(n-2:n, j-1:j+1, k-1:k+1, :));
            x_cur = x;
            x_cur(n, j, k, l) = x(n, j, k, l) + dxmag;
            F1 = F(x_cur(n-2:n, j-1:j+1, k-1:k+1, :));
            g(n, j, k, l) = (F1-F0)/dxmag;
        end
    end
end
for i = 2:n-1
    for k = 2:n-1
        for l = 1:6
            F0 = F(x(i-1:i+1, 1:3, k-1:k+1, :));
            x_cur = x;
            x_cur(i, 1, k, l) = x(i, 1, k, l) + dxmag;
            F1 = F(x_cur(i-1:i+1, 1:3, k-1:k+1, :));
            g(i, 1, k, l) = (F1-F0)/dxmag;

            F0 = F(x(i-1:i+1, n-2:n, k-1:k+1, :));
            x_cur = x;
            x_cur(i, n, k, l) = x(i, n, k, l) + dxmag;
            F1 = F(x_cur(i-1:i+1, n-2:n, k-1:k+1, :));
            g(i, n, k, l) = (F1-F0)/dxmag;
        end
    end
end
for i = 2:n-1
    for j = 2:n-1
        for l = 1:6
            F0 = F(x(i-1:i+1, j-1:j+1, 1:3, :));
            x_cur = x;
            x_cur(i, j, 1, l) = x(i, j, 1, l) + dxmag;
            F1 = F(x_cur(i-1:i+1, j-1:j+1, 1:3, :));
            g(i, j, 1, l) = (F1-F0)/dxmag;

            F0 = F(x(i-1:i+1, j-1:j+1, n-2:n, :));
            x_cur = x;
            x_cur(i, j, n, l) = x(i, j, n, l) + dxmag;
            F1 = F(x_cur(i-1:i+1, j-1:j+1, n-2:n, :));
            g(i, j, n, l) = (F1-F0)/dxmag;
        end
    end
end
for i = 2:n-1
    for l = 1:6
        F0 = F(x(i-1:i+1, 1:3, 1:3, :));
        x_cur = x;
        x_cur(i, 1, 1, l) = x(i, 1, 1, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, 1:3, 1:3, :));
        g(i, 1, 1, l) = (F1-F0)/dxmag;

        F0 = F(x(i-1:i+1, n-2:n, 1:3, :));
        x_cur = x;
        x_cur(i, n, 1, l) = x(i, n, 1, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, n-2:n, 1:3, :));
        g(i, n, 1, l) = (F1-F0)/dxmag;

        F0 = F(x(i-1:i+1, 1:3, n-2:n, :));
        x_cur = x;
        x_cur(i, 1, n, l) = x(i, 1, n, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, 1:3, n-2:n, :));
        g(i, 1, n, l) = (F1-F0)/dxmag;

        F0 = F(x(i-1:i+1, n-2:n, n-2:n, :));
        x_cur = x;
        x_cur(i, n, n, l) = x(i, n, n, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, n-2:n, n-2:n, :));
        g(i, n, n, l) = (F1-F0)/dxmag;
    end
end

for j = 2:n-1
    for l = 1:6
        F0 = F(x(1:3, j-1:j+1, 1:3, :));
        x_cur = x;
        x_cur(1, j, 1, l) = x(1, j, 1, l) + dxmag;
        F1 = F(x_cur(1:3, j-1:j+1, 1:3, :));
        g(1, j, 1, l) = (F1-F0)/dxmag;

        F0 = F(x(n-2:n, j-1:j+1, 1:3, :));
        x_cur = x;
        x_cur(n, j, 1, l) = x(n, j, 1, l) + dxmag;
        F1 = F(x_cur(n-2:n, j-1:j+1, 1:3, :));
        g(n, j, 1, l) = (F1-F0)/dxmag;

        F0 = F(x(1:3, j-1:j+1, n-2:n, :));
        x_cur = x;
        x_cur(1, j, n, l) = x(1, j, n, l) + dxmag;
        F1 = F(x_cur(1:3, j-1:j+1, n-2:n, :));
        g(1, j, n, l) = (F1-F0)/dxmag;

        F0 = F(x(n-2:n, j-1:j+1, n-2:n, :));
        x_cur = x;
        x_cur(n, j, n, l) = x(n, j, n, l) + dxmag;
        F1 = F(x_cur(n-2:n, j-1:j+1, n-2:n, :));
        g(n, j, n, l) = (F1-F0)/dxmag;
    end
end
for i = 2:n-1
    for l = 1:6
        F0 = F(x(i-1:i+1, 1:3, 1:3, :));
        x_cur = x;
        x_cur(i, 1, 1, l) = x(i, 1, 1, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, 1:3, 1:3, :));
        g(i, 1, 1, l) = (F1-F0)/dxmag;

        F0 = F(x(i-1:i+1, n-2:n, 1:3, :));
        x_cur = x;
        x_cur(i, n, 1, l) = x(i, n, 1, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, n-2:n, 1:3, :));
        g(i, n, 1, l) = (F1-F0)/dxmag;

        F0 = F(x(i-1:i+1, 1:3, n-2:n, :));
        x_cur = x;
        x_cur(i, 1, n, l) = x(i, 1, n, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, 1:3, n-2:n, :));
        g(i, 1, n, l) = (F1-F0)/dxmag;

        F0 = F(x(i-1:i+1, n-2:n, n-2:n, :));
        x_cur = x;
        x_cur(i, n, n, l) = x(i, n, n, l) + dxmag;
        F1 = F(x_cur(i-1:i+1, n-2:n, n-2:n, :));
        g(i, n, n, l) = (F1-F0)/dxmag;
    end
end
for k = 2:n-1
    for l = 1:6
        F0 = F(x(1:3, 1:3, k-1:k+1, :));
        x_cur = x;
        x_cur(1, 1, k, l) = x(1, 1, k, l) + dxmag;
        F1 = F(x_cur(1:3, 1:3, k-1:k+1, :));
        g(1, 1, k, l) = (F1-F0)/dxmag;

        F0 = F(x(n-2:n, 1:3, k-1:k+1, :));
        x_cur = x;
        x_cur(n, 1, k, l) = x(n, 1, k, l) + dxmag;
        F1 = F(x_cur(n-2:n, 1:3, k-1:k+1, :));
        g(n, 1, k, l) = (F1-F0)/dxmag;

        F0 = F(x(1:3, n-2:n, k-1:k+1, :));
        x_cur = x;
        x_cur(1, n, k, l) = x(1, n, k, l) + dxmag;
        F1 = F(x_cur(1:3, n-2:n, k-1:k+1, :));
        g(1, n, k, l) = (F1-F0)/dxmag;

        F0 = F(x(n-2:n, n-2:n, k-1:k+1, :));
        x_cur = x;
        x_cur(n, n, k, l) = x(n, n, k, l) + dxmag;
        F1 = F(x_cur(n-2:n, n-2:n, k-1:k+1, :));
        g(n, n, k, l) = (F1-F0)/dxmag;
    end
end
for l = 1:6
    F0 = F(x(1:3, 1:3, 1:3, :));
    x_cur = x;
    x_cur(1, 1, 1, l) = x(1, 1, 1, l) + dxmag;
    F1 = F(x_cur(1:3, 1:3, 1:3, :));
    g(1, 1, 1, l) = (F1-F0)/dxmag;

    F0 = F(x(n-2:n, 1:3, 1:3, :));
    x_cur = x;
    x_cur(n, 1, 1, l) = x(n, 1, 1, l) + dxmag;
    F1 = F(x_cur(n-2:n, 1:3, 1:3, :));
    g(n, 1, 1, l) = (F1-F0)/dxmag;

    F0 = F(x(1:3, n-2:n, 1:3, :));
    x_cur = x;
    x_cur(1, n, 1, l) = x(1, n, 1, l) + dxmag;
    F1 = F(x_cur(1:3, n-2:n, 1:3, :));
    g(1, n, 1, l) = (F1-F0)/dxmag;
    
    F0 = F(x(1:3, 1:3, n-2:n, :));
    x_cur = x;
    x_cur(1, 1, n, l) = x(1, 1, n, l) + dxmag;
    F1 = F(x_cur(1:3, 1:3, n-2:n, :));
    g(1, 1, n, l) = (F1-F0)/dxmag;
    
    F0 = F(x(n-2:n, n-2:n, 1:3, :));
    x_cur = x;
    x_cur(n, n, 1, l) = x(n, n, 1, l) + dxmag;
    F1 = F(x_cur(n-2:n, n-2:n, 1:3, :));
    g(n, n, 1, l) = (F1-F0)/dxmag;
    
    F0 = F(x(n-2:n, 1:3, n-2:n, :));
    x_cur = x;
    x_cur(n, 1, n, l) = x(n, 1, n, l) + dxmag;
    F1 = F(x_cur(n-2:n, 1:3, n-2:n, :));
    g(n, 1, n, l) = (F1-F0)/dxmag;
    
    F0 = F(x(1:3, n-2:n, n-2:n, :));
    x_cur = x;
    x_cur(1, n, n, l) = x(1, n, n, l) + dxmag;
    F1 = F(x_cur(1:3, n-2:n, n-2:n, :));
    g(1, n, n, l) = (F1-F0)/dxmag;
    
    F0 = F(x(n-2:n, n-2:n, n-2:n, :));
    x_cur = x;
    x_cur(n, n, n, l) = x(n, n, n, l) + dxmag;
    F1 = F(x_cur(n-2:n, n-2:n, n-2:n, :));
    g(n, n, n, l) = (F1-F0)/dxmag;
end
g = reshape(g, 1, []);
end