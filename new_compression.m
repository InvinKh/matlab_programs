%% 2D
clc 
clear

n=10; %Кол-во атомов в стержне
z=1;% Межатомное расстояние
K=10; %Коэф. жесткости
a10= 0.2; %Граничные условия (на сколько оттягиваем левый конец);
a1=0; %Гр.усл. левый конец закреплен

a = zeros(n, n, 2);
a(n, n, 1) = 3;
a(n, n, 2) = 3;
a
% a(1)=a1;
% a(10)=a10;


c11 = 1;
c12 = 0.5;
c44 = 0.5;


L = z*n;
dx = L/n;
dy = L/n;

uxx0 = diff(a(:, :, 1))/dx;
uyy0 = diff(a(:, :, 2))/dy;
uxy0 = 1/2 * (diff(a(:, :, 1))/dy + diff(a(:, :, 2))/dx);


F_cub_2d_arr = @(c11,c12,c44,uxx,uyy,uxy)   1/2 * c11*(uxx.^2+uyy.^2)   + c12*(uxx.*uyy) + 2*c44*(uxy.^2);


F_cub_2d = @(c11,c12,c44,uxx,uyy,uxy) sum(sum(F_cub_2d_arr(c11,c12,c44,uxx,uyy,uxy)));

a_arr_new = F_cub_2d_arr(c11,c12,c44,uxx0, uyy0, uxy0)

F = F_cub_2d(c11,c12,c44,uxx0, uyy0, uxy0)

%%
N = 10;
dx = L/N;
dy = L/N;
a_arr = zeros(N, N, 2);
a_arr_0 = reshape(a_arr, [1, 2*N*N]);

func_a_arr = @(a_arr, N) reshape(a_arr, [N, N, 2]);
func_uxx = @(a_arr, dx) diff(a_arr(:, :, 1), 1, 1)/dx;
func_uyy = @(a_arr, dy) (diff(a_arr(:, :, 2), 1, 2)/dy);
func_uxy_x = @(a_arr, dx) diff(a_arr(:, 1:N-1, 2), 1, 1)/dx;
func_uxy_y = @(a_arr, dy) diff(a_arr(1:N-1, :, 1), 1, 2)/dy;
func_uxy = @(a_arr, dx, dy) 1/2 * (abs(func_uxy_x(a_arr, dx)) + abs(func_uxy_y(a_arr, dy)));

F_cub_2d = @(c11,c12,c44,uxx,uyy,uxy) sum(sum(1/2 * c11*(uxx.^2))) + sum(sum(1/2 * c11*(uyy.^2)))   + sum(sum(c12*(uxx.*uyy))) + sum(sum(2*c44*(uxy.^2)));

optfunc = @(a_arr) F_cub_2d(c11,c12,c44,func_uxx(func_a_arr(a_arr, N),dx), func_uyy(func_a_arr(a_arr, N), dy), func_uxy(func_a_arr(a_arr, N), dx, dy));

u_xx = func_uxx(a_arr, dx)
u_yy = func_uyy(a_arr, dy)

u_xy_x = func_uxy_x(a_arr, dx)
u_xy_y = func_uxy_y(a_arr, dy)

x = linspace(-4, 4, N);


Aeq = zeros(2*N*N,2*N*N);
for i=0:(N-1)
    Aeq(i*N+1,i*N+1) = 1;
    Aeq(N*N+i*N+1,N*N+i*N+1) = 1;
end

% Aeq(1,1) = 1; 
% Aeq(N*N+1,N*N+1) = 1;
% Aeq(N,N) = 1; 
% Aeq(N*N+N,N*N+N) = 1;
% % Aeq(N*(N-1)+1,N*(N-1)+1) = 1; 
% % Aeq(N*N+N*(N-1)+1,N*N+N*(N-1)+1) = 1;
% % Aeq(N*N,N*N) = 1; 
% % Aeq(2*N*N,2*N*N) = 1;



beq = zeros(2*N*N, 1);
for i=0:(N-1)
    beq(i*N+1, 1) = 0;
    beq(N*N+i*N+1, 1) = x(i+1);
end

% beq = zeros(2*N*N, 1); 
% beq(1, 1) = -7;
% beq(N, 1) = 7; 
% beq(N*(N-1)+1, 1) = 0;
% beq(N*N, 1) = 0;
% %Это не нужно, но я для себя написал Г.У. Можно закоментить.
% beq(N*N+N, 1) = 0; 
% beq(N*N+N*(N-1)+1, 1) = 0;
% beq(2*N*N, 1) = 0;
% Aeq

%x - смещение по ox вправо и влево. y - смещение по oy вверх и вниз
% for i=1:N
%     beq()
% end

lb = [];
ub = [];
nonlcon = [];

% options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000,'Algorithm', 'sqp','Display','iter', 'MaxIterations', 3000);
options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000, 'Algorithm', 'sqp', 'Display','iter', 'MaxIterations', 3000);
% 
tic
[a_arr1, fmin]= fmincon(optfunc,a_arr_0,[],[],Aeq,beq, [],[],[],options);
toc

a_arr1 = func_a_arr(a_arr1, N);
Cord=[]
for i=1:N
for j=1:N
Cord(i,j,1)=i+a_arr1(i,j,1);
Cord(i,j,2)=j+a_arr1(i,j,2);
end
end

figure
hold on
for i=1:N
for j=1:N
plot(Cord(i,j,1), Cord(i,j,2), 'ro');
end
end
% figure
% hold on
% grid on
% title("N = " + N);
% plot(a_arr1(:, :, 1), a_arr1(:, :, 2), "DisplayName", ['N = ' N]);
%%
N = 10;

func = optifunc(a_arr);
% f = fmincon(optifunc,a_arr_0,[],[],Aeq,beq, [],[],[],options);
%%
function f = optifunc(input_arr)
    N = (numel(input_arr)/2)^(0.5);
    a_arr = input_arr(1:end);
    a_arr = reshape(a_arr, [N, N, 2]);
    [uxx, uyy] = gradient(a_arr(:, :, 1));

    f = sum(sum(sum(a_arr)));
end
