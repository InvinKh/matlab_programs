%% Interaction 3D
clc 
clear

N=3; %Кол-во атомов в стержне

elastic_arr = zeros(N, N, N, 3);
% dom_arr = rand(N, N, N, 3);
dom_arr = zeros(N, N, N, 3);
work_arr = zeros(N, N, N, 6);
work_arr(:, :, :, 1:3) = elastic_arr;
work_arr(:, :, :, 4:6) = dom_arr;

opt_arr = reshape(work_arr, [1, 6*N*N*N]);

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


% sum(A, 'all'); fix
energy_domen_3d = @(domen_arr) sum(sum(sum( ...
    domen_arr_3d(a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, ...
    a1123, domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), ...
    domen_arr(:, :, :, 3)))));

energy_elasticity_3d = @(a_arr) sum(sum(sum(elasticity_arr_3d(c11,c12,c44,...
    func_u_3d(func_arr_3d(a_arr, N), 1, 1, dx, dx), ...
    func_u_3d(func_arr_3d(a_arr, N), 2, 2, dy, dy), ...
    func_u_3d(func_arr_3d(a_arr, N), 3, 3, dz, dz), ...
    func_u_3d(func_arr_3d(a_arr, N), 1, 2, dx, dy), ...
    func_u_3d(func_arr_3d(a_arr, N), 2, 3, dy, dz), ...
    func_u_3d(func_arr_3d(a_arr, N), 1, 3, dx, dz)...
    ))));

energy_interaction_3d = @(a_arr, domen_arr) sum(sum(sum(interaction_arr_3d(q11, q12, q44, ...
    func_u_3d(a_arr, 1, 1, dx, dx), ...
    func_u_3d(a_arr, 2, 2, dy, dy), ...
    func_u_3d(a_arr, 3, 3, dz, dz), ...
    func_u_3d(a_arr, 1, 2, dx, dy), ...
    func_u_3d(a_arr, 2, 3, dy, dz), ...
    func_u_3d(a_arr, 1, 3, dx, dz), ...
    domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), domen_arr(:, :, :, 3) ...
    ))));

energy_domen_gradient_3d = @(domen_arr) sum(sum(sum(domen_gradient_arr_3d(G11, G14, G44, ...
    domen_diff_3d(domen_arr(:, :, :, 1), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 3, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 1), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 3, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 1), 3, dx) ...
    ))));

optfunc = @(a_arr) energy_elasticity_3d(func_arr_3d(a_arr(1:numel(a_arr)/2), N)) ... 
    + energy_domen_3d(func_arr_3d(a_arr(numel(a_arr)/2+1:end), N)) ...
    + energy_interaction_3d(func_arr_3d(a_arr(1:numel(a_arr)/2), N), ...
                            func_arr_3d(a_arr(numel(a_arr)/2+1:end), N)) ...
    + energy_domen_gradient_3d(func_arr_3d(a_arr(numel(a_arr)/2+1:end), N));

% Set the boundary conditions
Aeq = zeros(6*N*N*N,6*N*N*N);
beq = zeros(6*N*N*N, 1);

% % Corners of the lower face (only elastic)
% % Points on axe x
% Aeq(1, 1) = 1; Aeq(N, N) = 1;
% Aeq(N*(N-1)+1, N*(N-1)+1) = 1; Aeq(N*N, N*N) = 1;
% % Points on axe y
% Aeq(N*N*N+1, N*N*N+1) = 1; Aeq(N*N*N+N, N*N*N+N) = 1;
% Aeq(N*N*N+N*(N-1)+1, N*N*N+N*(N-1)+1) = 1; Aeq(N*N*N+N*N, N*N*N+N*N) = 1;
% % Points on axe z
% Aeq(2*N*N*N+1, 2*N*N*N+1) = 1; Aeq(2*N*N*N+N, 2*N*N*N+N) = 1;
% Aeq(2*N*N*N+N*(N-1)+1, 2*N*N*N+N*(N-1)+1) = 1; Aeq(2*N*N*N+N*N, 2*N*N*N+N*N) = 1;
% 
% beq(1, 1) = -0.5; beq(N, 1) = 0.5; beq(N*(N-1)+1, 1) = -0.5; beq(N*N) = 0.5; 
% beq(N*N*N+1, 1) = -0.5; beq(N*N*N+N, 1) = -0.5; beq(N*N*N+N*(N-1)+1, 1) = 0.5; beq(N*N*N+N*N) = 0.5; 
% beq(2*N*N*N+1, 1) = 0; beq(2*N*N*N+N, 1) = 0; beq(2*N*N*N+N*(N-1)+1, 1) = 0; beq(2*N*N*N+N*N) = 0;

% All lower face
% Corners on axe x
Aeq(1, 1) = 1; Aeq(N, N) = 1;
Aeq(N*(N-1)+1, N*(N-1)+1) = 1; Aeq(N*N, N*N) = 1;
% Corners on axe y
Aeq(N*N*N+1, N*N*N+1) = 1; Aeq(N*N*N+N, N*N*N+N) = 1;
Aeq(N*N*N+N*(N-1)+1, N*N*N+N*(N-1)+1) = 1; Aeq(N*N*N+N*N, N*N*N+N*N) = 1;
% Points on axe z (all lower face)
for i = 1:N*N
    Aeq(2*N*N*N+i, 2*N*N*N+i) = 1;
end

beq(1, 1) = -0.5; beq(N, 1) = 0.5; beq(N*(N-1)+1, 1) = -0.5; beq(N*N) = 0.5; 
beq(N*N*N+1, 1) = -0.5; beq(N*N*N+N, 1) = -0.5; beq(N*N*N+N*(N-1)+1, 1) = 0.5; beq(N*N*N+N*N) = 0.5; 
for i = 1:N*N
    beq(2*N*N*N+i, 1) = 0;
end
 




options = optimoptions(@fmincon, 'Algorithm', 'sqp','Display','iter', 'MaxFunctionEvaluations', 10^7, 'MaxIterations', 10^3, 'UseParallel', true);
% options = optimoptions(@fmincon,'Algorithm', 'sqp','Display','iter', 'MaxIterations', 1000);
% options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations', 10, 'Algorithm', 'interior-point', 'UseParallel', true);
% options.SubproblemAlgorithm = "cg";
% options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations', 100, 'Algorithm', 'trust-region-reflective');

tic
[a_arr1,fval,exitflag,output,lambda,grad,hessian] = fmincon(optfunc,work_arr,[],[],Aeq,beq, [],[],[],options);
toc

%%
figure;
surf(hessian)
%%
figure;
plot(grad)

%% Time/accuracy
for N=4:7
for X=50:50:300

elastic_arr = zeros(N, N, N, 3);
% dom_arr = rand(N, N, N, 3);
dom_arr = zeros(N, N, N, 3);
work_arr = zeros(N, N, N, 6);
work_arr(:, :, :, 1:3) = elastic_arr;
work_arr(:, :, :, 4:6) = dom_arr;

opt_arr = reshape(work_arr, [1, 6*N*N*N]);

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

energy_domen_3d = @(domen_arr) sum(sum(sum( ...
    domen_arr_3d(a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, ...
    a1123, domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), ...
    domen_arr(:, :, :, 3)))));

energy_elasticity_3d = @(a_arr) sum(sum(sum(elasticity_arr_3d(c11,c12,c44,...
    func_u_3d(func_arr_3d(a_arr, N), 1, 1, dx, dx), ...
    func_u_3d(func_arr_3d(a_arr, N), 2, 2, dy, dy), ...
    func_u_3d(func_arr_3d(a_arr, N), 3, 3, dz, dz), ...
    func_u_3d(func_arr_3d(a_arr, N), 1, 2, dx, dy), ...
    func_u_3d(func_arr_3d(a_arr, N), 2, 3, dy, dz), ...
    func_u_3d(func_arr_3d(a_arr, N), 1, 3, dx, dz)...
    ))));

energy_interaction_3d = @(a_arr, domen_arr) sum(sum(sum(interaction_arr_3d(q11, q12, q44, ...
    func_u_3d(a_arr, 1, 1, dx, dx), ...
    func_u_3d(a_arr, 2, 2, dy, dy), ...
    func_u_3d(a_arr, 3, 3, dz, dz), ...
    func_u_3d(a_arr, 1, 2, dx, dy), ...
    func_u_3d(a_arr, 2, 3, dy, dz), ...
    func_u_3d(a_arr, 1, 3, dx, dz), ...
    domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), domen_arr(:, :, :, 3) ...
    ))));

energy_domen_gradient_3d = @(domen_arr) sum(sum(sum(domen_gradient_arr_3d(G11, G14, G44, ...
    domen_diff_3d(domen_arr(:, :, :, 1), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 3, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 1), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 2), 3, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 2, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 3), 1, dx), ...
    domen_diff_3d(domen_arr(:, :, :, 1), 3, dx) ...
    ))));

optfunc = @(a_arr) energy_elasticity_3d(func_arr_3d(a_arr(1:numel(a_arr)/2), N)) ... 
    + energy_domen_3d(func_arr_3d(a_arr(numel(a_arr)/2+1:end), N)) ...
    + energy_interaction_3d(func_arr_3d(a_arr(1:numel(a_arr)/2), N), ...
                            func_arr_3d(a_arr(numel(a_arr)/2+1:end), N)) ...
    + energy_domen_gradient_3d(func_arr_3d(a_arr(numel(a_arr)/2+1:end), N));

% Set the boundary conditions
Aeq = zeros(6*N*N*N,6*N*N*N);
beq = zeros(6*N*N*N, 1);

% % Corners of the lower face (only elastic)
% % Points on axe x
% Aeq(1, 1) = 1; Aeq(N, N) = 1;
% Aeq(N*(N-1)+1, N*(N-1)+1) = 1; Aeq(N*N, N*N) = 1;
% % Points on axe y
% Aeq(N*N*N+1, N*N*N+1) = 1; Aeq(N*N*N+N, N*N*N+N) = 1;
% Aeq(N*N*N+N*(N-1)+1, N*N*N+N*(N-1)+1) = 1; Aeq(N*N*N+N*N, N*N*N+N*N) = 1;
% % Points on axe z
% Aeq(2*N*N*N+1, 2*N*N*N+1) = 1; Aeq(2*N*N*N+N, 2*N*N*N+N) = 1;
% Aeq(2*N*N*N+N*(N-1)+1, 2*N*N*N+N*(N-1)+1) = 1; Aeq(2*N*N*N+N*N, 2*N*N*N+N*N) = 1;
% 
% beq(1, 1) = -0.5; beq(N, 1) = 0.5; beq(N*(N-1)+1, 1) = -0.5; beq(N*N) = 0.5; 
% beq(N*N*N+1, 1) = -0.5; beq(N*N*N+N, 1) = -0.5; beq(N*N*N+N*(N-1)+1, 1) = 0.5; beq(N*N*N+N*N) = 0.5; 
% beq(2*N*N*N+1, 1) = 0; beq(2*N*N*N+N, 1) = 0; beq(2*N*N*N+N*(N-1)+1, 1) = 0; beq(2*N*N*N+N*N) = 0;

% All lower face
% Corners on axe x
Aeq(1, 1) = 1; Aeq(N, N) = 1;
Aeq(N*(N-1)+1, N*(N-1)+1) = 1; Aeq(N*N, N*N) = 1;
% Corners on axe y
Aeq(N*N*N+1, N*N*N+1) = 1; Aeq(N*N*N+N, N*N*N+N) = 1;
Aeq(N*N*N+N*(N-1)+1, N*N*N+N*(N-1)+1) = 1; Aeq(N*N*N+N*N, N*N*N+N*N) = 1;
% Points on axe z (all lower face)
for i = 1:N*N
    Aeq(2*N*N*N+i, 2*N*N*N+i) = 1;
end

beq(1, 1) = -0.5; beq(N, 1) = 0.5; beq(N*(N-1)+1, 1) = -0.5; beq(N*N) = 0.5; 
beq(N*N*N+1, 1) = -0.5; beq(N*N*N+N, 1) = -0.5; beq(N*N*N+N*(N-1)+1, 1) = 0.5; beq(N*N*N+N*N) = 0.5; 
for i = 1:N*N
    beq(2*N*N*N+i, 1) = 0;
end
 




options = optimoptions(@fmincon, 'Algorithm', 'sqp','Display','iter', 'MaxFunctionEvaluations', 10^7, 'MaxIterations', 10^3, 'UseParallel', true);
% options = optimoptions(@fmincon,'Algorithm', 'sqp','Display','iter', 'MaxIterations', 1000);
% options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations',X, 'Algorithm', 'interior-point', 'UseParallel', true);
options.HessianFcn = [];
% 
tStart = tic;

[a_arr1,fval,exitflag,output,lambda,grad,hessian] = fmincon(optfunc,work_arr,[],[],Aeq,beq, [],[],[],options);
toc;
a = duration([0, 0, toc(tStart)]); 
number_point = string(N);
number_iteration = string(X);
name_direction = 'C:/Users/User/Labs/Matlab/resilience/new/';
name_file = name_direction + 'data' + number_point + '_' + number_iteration + '.mat';
save(name_file, 'a', 'a_arr1')
end
end
%% Draw 3D
value_arr_x = zeros(N, N, N);
value_arr_y = zeros(N, N, N);
value_arr_z = zeros(N, N, N);
value_arr_domen_x = zeros(N, N, N);
value_arr_domen_y = zeros(N, N, N);
value_arr_domen_z = zeros(N, N, N);
for i = 1:N
    value_arr_x(i, :, :) = a_arr1(i, :, :, 1) + i-1;
    value_arr_y(:, i, :) = a_arr1(:, i, :, 2) + i-1;
    value_arr_z(:, :, i) = a_arr1(:, :, i, 3) + i-1;
    value_arr_domen_x(i, :, :) = a_arr1(i, :, :, 4);
    value_arr_domen_y(:, i, :) = a_arr1(:, i, :, 5);
    value_arr_domen_z(:, :, i) = a_arr1(:, i, :, 6);
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
                value_arr_domen_x(i, j, k), value_arr_domen_y(i, j, k), value_arr_domen_z(i, j, k), 0.2);
%             plot3(value_arr_x(i, j, 1), value_arr_y(i, j, 1), value_arr_z(i, j, 1), 'o');
        end
    end
end
hold off;

%%
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

function array = domen_arr_3d(a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, a1123, p1, p2, p3)
    array = a1*(p1.^2 + p2.^2 + p3.^2) ...
          + a11*(p1.^4 + p2.^4 + p3.^4) ...
          + a12*(p1.*p2 + p2.*p3 + p1.*p3) ...
          + a111*(p1.^6 + p2.^6 + p3.^6) ...
          + a112*(p1.^4.*(p2.^2 + p3.^2) + p2.^4.*(p1.^2 + p3.^2) + p3.^4.*(p1.^2 + p2.^2)) ...
          + a123*(p1.^2.*p2.^2.*p3.^2) ...
          + a1111*(p1.^8 + p2.*8 + p3.*8) ...
          + a1112*(p1.^6.*(p2.^2 + p3.^2) + p2.^6.*(p1.^2 + p3.^2) + p3.^6.*(p1.^2 + p2.^2)) ...
          + a1122*(p1.^4.*p2.^4 + p2.^4.*p3.^4 + p1.^4.*p3.^4) ...
          + a1123*(p1.^4.*p2.^2.*p3.^2 + p1.^2.*p2.^4.*p3.^2 + p1.^2.*p2.^2.*p3.^4);
end

function array = elasticity_arr_3d(c11,c12,c44,u11,u22,u33, u12, u23, u13)
    array = 1/2 * c11*(u11.^2 + u22.^2 + u33.^2) ...
            + c12*(u11.*u22 + u22.*u33 + u11.*u33) ...
            + 2*c44*(u12.^2 + u23.^2 + u13.^2);
end

function array = interaction_arr_3d(q11, q12, q44, u11, u22, u33, u12, u23, u13, p1, p2, p3)
    array = -q11*(u11.*p1.^2 + u22.*p2.^2 + u33.*p3.^2)...
            -q12*(u11.*(p2.^2+p3.^2) + u22.*(p1.^2+p3.^2) + u33.*(p1.^2+p2.^2))...
            -2*q44*(u12.*p1.*p2 + u23.*p2.*p3 + u13.*p1.*p3);
end

function array = domen_gradient_arr_3d(G11, G14, G44, p11, p22, p33, p12, p21, p23, p32, p31, p13)
    array = 1/2*G11*(p11.^2 + p22.^2 + p33.^2) ...
            + G14*(p11.*p22 + p22.*p33 + p11.*p33) ...
            + G44*(p12.^2 + p21.^2 + p23.^2 + p32.^2 + p31.^2 + p13.^2);
end

% function [fval, grad] = energy_system(arr, a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, a1123)
% 
%     energy_domen_3d = @(domen_arr) sum(sum(sum( ...
%         domen_arr_3d(a1, a11, a12, a111, a112, a123, a1111, a1112, a1122, ...
%         a1123, domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), ...
%         domen_arr(:, :, :, 3)))));
% 
%     energy_elasticity_3d = @(a_arr) sum(sum(sum(elasticity_arr_3d(c11,c12,c44,...
%         func_u_3d(func_arr_3d(a_arr, N), 1, 1, dx, dx), ...
%         func_u_3d(func_arr_3d(a_arr, N), 2, 2, dy, dy), ...
%         func_u_3d(func_arr_3d(a_arr, N), 3, 3, dz, dz), ...
%         func_u_3d(func_arr_3d(a_arr, N), 1, 2, dx, dy), ...
%         func_u_3d(func_arr_3d(a_arr, N), 2, 3, dy, dz), ...
%         func_u_3d(func_arr_3d(a_arr, N), 1, 3, dx, dz)...
%         ))));
% 
%     energy_interaction_3d = @(a_arr, domen_arr) sum(sum(sum(interaction_arr_3d(q11, q12, q44, ...
%         func_u_3d(a_arr, 1, 1, dx, dx), ...
%         func_u_3d(a_arr, 2, 2, dy, dy), ...
%         func_u_3d(a_arr, 3, 3, dz, dz), ...
%         func_u_3d(a_arr, 1, 2, dx, dy), ...
%         func_u_3d(a_arr, 2, 3, dy, dz), ...
%         func_u_3d(a_arr, 1, 3, dx, dz), ...
%         domen_arr(:, :, :, 1), domen_arr(:, :, :, 2), domen_arr(:, :, :, 3) ...
%         ))));
% 
%     energy_domen_gradient_3d = @(domen_arr) sum(sum(sum(domen_gradient_arr_3d(G11, G14, G44, ...
%         domen_diff_3d(domen_arr(:, :, :, 1), 1, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 2), 2, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 3), 3, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 1), 2, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 2), 1, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 2), 3, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 3), 2, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 3), 1, dx), ...
%         domen_diff_3d(domen_arr(:, :, :, 1), 3, dx) ...
%         ))));
% 
% end

% Деление двух чисел друг на друга не даёт научной публикации