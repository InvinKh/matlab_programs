%вычисление координат атомов в стержне при известном давлении на правый
%конец при закрепленном левом конце
clear
clc

n=10; % Число "атомов" в стержне
K=7; % модуль жесткости стержня

x=linspace(0,n-1, n);
u=zeros(1,n);
s=zeros(1,n);

p=3; % Давление на правый конец

s=s-p; % Формула 2.7 и 2.9 для 1D случая превращаются в однородность давления по всему стержню
u=1/(K)*s; %Формула 4.7 (в знаменателе 3 меняется на 1 т.к. размерность 1D)

x1=x.*sqrt(1+2*u); % Формула в начале стр.11

figure
title(['p=', num2str(p), '     n=', num2str(n), '      K=', num2str(K)])
hold on
plot(x, zeros(1,n)+1, 'ro')
plot(x1, zeros(1,n)+1.5, 'bo')
ylim([0 2])
legend('Rod without pressure','Rod with pressure')
%% горцы вычисление энергии
clc 
clear

a=rand(1,10); 
% a=[0.1 0.15 0.05 0.07 0.11 0.12 0.12 0.05 0.1 0.14];

K=4;

for i=1:numel(a)-1
    u(i)=a(i+1)-a(i);
end

f= K*u.*u;

F=sum(f)

%% задача про растяжение/сжатие стержня на определенную величину
clc 
clear

n=10; %Кол-во атомов в стержне
z=1;% Межатомное расстояние
K=10; %Коэф. жесткости
a10= 11; %Граничные условия (на сколько оттягиваем левый конец);
a1=-11; %Гр.усл. левый конец закреплен

a=zeros(1,n);
a(1)=a1;
a(10)=a10;

F= @(y)K*((y(2)-a(1))^2+(y(3)-y(2))^2+(y(4)-y(3))^2+(y(5)-y(4))^2+(y(6)-y(5))^2+(y(7)-y(6))^2+(y(8)-y(7))^2+(y(9)-y(8))^2+(a(10)-y(9))^2);
% F= @(y)K*((y(2)-y(1))^2+(y(3)-y(2))^2+(y(4)-y(3))^2+(y(5)-y(4))^2+(y(6)-y(5))^2+(y(7)-y(6))^2+(y(8)-y(7))^2+(y(9)-y(8))^2+(y(10)-y(9))^2);


[y, fmin]= fminunc(F,a);

for i=1:n
    x(i)=(i-1)*z+y(i);
    x1(i)=(i-1)*z;
end

figure
hold on
title(['a10=', num2str(a10), '     n=', num2str(n-1), '      z=', num2str(z), '        F=', num2str(fmin)])
plot(x1, zeros(1,n)+1.5, 'ro')
plot(x, zeros(1,n)+1, 'bo')
ylim([0 2])
legend('Rod without pressure','Rod under the pressure')

%% задача как в первой секции, только через горцев (не работает)
clc 
clear




n=10; %Кол-во атомов в стержне
u=zeros(1,n-1);
a=zeros(1, n);
z=1;% Межатомное расстояние
K=7; %Коэф. жесткости
p=3 % Давление на правый конец
a1=0; %Гр.усл. левый конец закреплен

u(n-1)= -p/K;


F= @(y)abs((y(2)-a(1))+(y(3)-y(2))+(y(4)-y(3))+(y(5)-y(4))+(y(6)-y(5))+(y(7)-y(6))+(y(8)-y(7))+(y(9)-y(8))+(u(n-1)));

%abs((y(2)-a(1))+(y(3)-y(2))+(y(4)-y(3))+(y(5)-y(4))+(y(6)-y(5))+(y(7)-y(6))+(y(8)-y(7))+(y(9)-y(8))+(u(n-1)))*
[y, fmin]= fminunc(F,a);

for i=1:n
    x(i)=(i-1)*z+y(i);
    x1(i)=(i-1)*z;
end

figure
hold on
title(['     n=', num2str(n-1), '      z=', num2str(z),'     K=', num2str(K)])
plot(x1, zeros(1,n)+1, 'ro')
plot(x, zeros(1,n)+1.5, 'bo')
ylim([0 2])
legend('Rod without pressure','Rod with pressure')



%% 1D
clc 
clear

n=10; %Кол-во атомов в стержне
z=1;% Межатомное расстояние
K=10; %Коэф. жесткости
a10= 0.2; %Граничные условия (на сколько оттягиваем левый конец);
a1=0; %Гр.усл. левый конец закреплен

a=zeros(1,n);
a(1)=a1;
a(10)=a10;


c11 = 1;
c12 = 0.5;
c44 = 0.5;

ax = linspace(a1,a10,n)

L = 10;
dx = L/n;

uxx0 = diff(ax)/dx


F_cub_2d_arr = @(c11,c12,c44,uxx,uyy,uxy)   1/2 * c11*(uxx.^2+uyy.^2)   + c12*(uxx.*uyy) + 2*c44*(uxy.^2);


F_cub_2d = @(c11,c12,c44,uxx,uyy,uxy) sum(sum(F_cub_2d_arr(c11,c12,c44,uxx,uyy,uxy)));

F_cub_2d_arr(c11,c12,c44,uxx0, uxx0*0, uxx0*0)

F_cub_2d(c11,c12,c44,uxx0, uxx0*0, uxx0*0)

%% Minimize it (1D)
N = 50;
dx = L/N;
a_arr_0 = linspace(0,1,N)
func_uxx = @(a_arr, dx) diff(a_arr)/dx;

optfunc = @(a_arr) F_cub_2d(c11,c12,c44,func_uxx(a_arr,dx), func_uxx(a_arr,dx)*0, func_uxx(a_arr,dx)*0);
checkfunc = @() fmincon(optfunc,a_arr_0,[],[],Aeq,beq, lb,ub,nonlcon,options);
%[a_arr1, fmin]= fminunc(optfunc,a_arr_0)

Aeq = zeros(N,N); Aeq(1,1) = 1; Aeq(N,N) = 1;
beq = zeros(N,1); beq(N) = 3;

lb = [];
ub = [];
nonlcon = [];

% options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000,'Algorithm', 'sqp','Display','iter', 'MaxIterations', 3000);
options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000, 'Algorithm', 'sqp', 'Display','iter', 'MaxIterations', 3000);

tic
[a_arr1, fmin]= fmincon(optfunc,a_arr_0,[],[],Aeq,beq, lb,ub,nonlcon,options)
toc

figure
hold on
grid on
title("N = " + N);
plot(a_arr1, "DisplayName", ['N = ' N]);
%% time/number dependence
numIt = 10;
N = zeros(1,numIt);
t = zeros(1,numIt);
for i=1:numIt
    N(i) = i*100;
    dx = L/N(i);
    a_arr_0 = linspace(0,1,N(i));
    func_uxx = @(a_arr, dx) diff(a_arr)/dx;
    optfunc = @(a_arr) F_cub_2d(c11,c12,c44,func_uxx(a_arr,dx), func_uxx(a_arr,dx)*0, func_uxx(a_arr,dx)*0);
    Aeq = zeros(N(i),N(i)); Aeq(1,1) = 1; Aeq(N(i),N(i)) = 1;
    beq = zeros(N(i),1); beq(N(i)) = 3;
    lb = [];
    ub = [];
    nonlcon = [];
    checkfunc = @() fmincon(optfunc,a_arr_0,[],[],Aeq,beq, lb,ub,nonlcon,options);
    options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000,'Algorithm', 'sqp', 'MaxIterations', 3000);
    t(i) = timeit(checkfunc);
end
figure
hold on
grid on
title("dependence of time on quantity");
xlabel("N", "Fontsize",18)
ylabel("t, second","Fontsize",18)

plot(N(1:end), t(1:end), "DisplayName", ['N = ' N]);
savefig(['C:\Users\User\Labs\Matlab\resilience\temedep.fig'])
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

%% !!!NEW!!! Minimize 
N = 10;
dx = L/N;
dy = L/N;
a_arr = zeros(N, N, 2);
a_arr_0 = reshape(a_arr, [1, 2*N*N]);

func_arr = @(a_arr, N) reshape(a_arr, [N, N, 2]);
func_uxx = @(a_arr, dx) diff(a_arr(:, 1:N-1, 1), 1, 1)/dx;
func_uyy = @(a_arr, dy) (diff(a_arr(1:N-1, :, 2), 1, 2)/dy);
func_uxy_x = @(a_arr, dx) diff(a_arr(:, 1:N-1, 2), 1, 1)/dx;
func_uxy_y = @(a_arr, dy) diff(a_arr(1:N-1, :, 1), 1, 2)/dy;
func_uxy = @(a_arr, dx, dy) 1/2 * (abs(func_uxy_x(a_arr, dx)) + abs(func_uxy_y(a_arr, dy)));

F_cub_2d = @(c11,c12,c44,uxx,uyy,uxy) sum(sum(1/2 * c11*(uxx.^2))) + sum(sum(1/2 * c11*(uyy.^2)))   + sum(sum(c12*(uxx.*uyy))) + sum(sum(2*c44*(uxy.^2)));

optfunc = @(a_arr) F_cub_2d(c11,c12,c44,func_uxx(func_arr(a_arr, N),dx), func_uyy(func_arr(a_arr, N), dy), func_uxy(func_arr(a_arr, N), dx, dy));


x = linspace(-1, 1, N);


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
options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000, 'Algorithm', 'sqp', 'Display','iter', 'MaxIterations', 1000);
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


%% Test

a_arr = randi([0 10], 4, 4)
l = 1;
nd_arr=zeros(numel(a_arr(:,1)), numel(a_arr(1,:)));
switch l
    case 1
        d_arr = diff(a_arr, 1, 2);
        for i=1:numel(a_arr(:,1))
            for j=1:numel(a_arr(1,:))-1
                nd_arr(i,j)=nd_arr(i,j) + d_arr(i,j)/2;
                nd_arr(i,j+1) = nd_arr(i,j+1) + d_arr(i,j)/2;
            end
        end
    case 2
        d_arr = diff(a_arr, 1, 1);
        for j=1:numel(a_arr(1,:))
            for i=1:numel(a_arr(:,1))-1
                nd_arr(i,j)=nd_arr(i,j) + d_arr(i,j)/2;
                nd_arr(i+1,j) = nd_arr(i+1,j) + d_arr(i,j)/2;
            end
        end
end
diff(a_arr, 1, 2)
nd_arr

%% test gradient
a_arr = randi([0 10], 4, 4)
[grad_xx, grad_xy] = gradient(a_arr)

%% Domen 1D
clc 
clear

N=100; %Кол-во атомов в стержне

L = 1;

v_arr = rand(N, 1);

A = -0.7; %A<0
B = 0.3; %B>0
D = 0.001;


dx = L/N;

F_domen_cub_1d_arr = @(A, B, v)   A*(v.^2) + B*(v.^4);

F_domen_cub_1d = @(A, B, v) sum(F_domen_cub_1d_arr(A, B, v));

F_domen_gradient_1d_arr = @(D, dx,  v) D*(diff(v)/dx).^2;

F_domen_gradient_1d =  @(D, dx, v) sum(F_domen_gradient_1d_arr(D, dx, v));

domen_optfunc = @(v) F_domen_cub_1d(A, B,v) + F_domen_gradient_1d(D, dx, v);

options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations', 200);

Aeq = zeros(N,1); Aeq(1,1) = 1; Aeq(N,N) = 1;
beq = zeros(N,1); beq(1) = 1.08; beq(N) = -1.08;
% 
% F = domen_optfunc(domen_arr)
% 
tic
[a_arr1, fmin]= fmincon(domen_optfunc,v_arr,[],[],Aeq,beq,[],[],[],options)
toc

xx = 1:1:N;

figure;
plot(xx, a_arr1);

%% Interaction 1D
clc 
clear

N=20; %Кол-во атомов в стержне

L = 1;
dx = L/N;

a=zeros(N, 1);

% v_arr = rand(N, 1);
v_arr = zeros(N, 1);

opt_arr = zeros(2*N, 1);
opt_arr(1:N, 1) = a;
opt_arr(N+1:end, 1) = v_arr;

c11 = 1;
c12 = 0.5;
c44 = 0.5;

A = -0.7; %A<0
B = 0.3; %B>0
D = 0.001;
E = 0.5;

uxx0 = zeros(N-1, 1);
uxx = @(a, dx) diff(a)/dx;


F_cub_1d_arr = @(c11,c12,c44,a, dx)   1/2*c11*(uxx(a, dx).^2);


F_cub_1d = @(c11,c12,c44,a,dx) sum(F_cub_1d_arr(c11,c12,c44,a,dx));
% 
% F_cub_1d_arr(c11,c12,c44,uxx0, uxx0*0, uxx0*0)
% 
% F_cub_1d(c11,c12,c44,uxx0, uxx0*0, uxx0*0)





F_domen_cub_1d_arr = @(A, B, v)   A*(v.^2) + B*(v.^4);

F_domen_cub_1d = @(A, B, v) sum(F_domen_cub_1d_arr(A, B, v));

F_domen_gradient_1d_arr = @(D, dx,  v) D*(diff(v)/dx).^2;

F_domen_gradient_1d =  @(D, dx, v) sum(F_domen_gradient_1d_arr(D, dx, v));

domen_optfunc = @(v) F_domen_cub_1d(A, B,v) + F_domen_gradient_1d(D, dx, v);

F_dependence_arr = @(E, v, a, dx) E*v.*uxx(a, dx);

F_dependence = @(E, v, a, dx) sum(F_dependence_arr(E, v, a, dx));

opt_func = @(arr) domen_optfunc(arr(N+1:end, 1)) + F_cub_1d(c11,c12,c44,arr(1:N, 1), dx) + F_dependence(E, arr(N+2:end, 1), arr(1:N, 1), dx);

options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations', 1300, 'OutputFcn', @domen_plot);

Aeq = zeros(2*N,2*N); Aeq(1,1) = 1; Aeq(N,N) = 1; Aeq(N+1,N+1) = 1; Aeq(2*N,2*N) = 1;
beq = zeros(2*N,1); beq(1) = 0; beq(N) = 0.3; beq(N+1) = 1.08; beq(2*N) = -1.08;
% 
% F = domen_optfunc(domen_arr)
% 
figure;
tic
[a_arr1, fmin]= fmincon(opt_func,opt_arr,[],[],Aeq,beq,[],[],[],options)
toc

xx = 1:1:N;
% 

%%
xx_u = 1:1:N-1;
figure;
plot(xx, a_arr1(1:N, 1));
% figure;
% plot(xx, a_arr1(N+1:end, 1));
% figure;
% plot(xx_u, uxx(a_arr1(1:N, 1), dx))

%% Domen 2D
clc 
clear

N=5; %Кол-во атомов в стержне
z=1;% Межатомное расстояние

b_arr = rand(N, N, 2);
% b_arr = zeros(N, N, 2);
% b_arr(1, 1, 1) = 0.008;
% b_arr(1, 1, 2) = 0;
% b_arr = -1 + 2*rand(N, N, 2);
b_arr_0 = reshape(b_arr, [1, 2*N*N]);


Ax = -0.7; %A<0
Ay = -0.7; %A<0
Bx = 0.3; %B>0
By = 0.3; %B>0


L = z*N;
dx = L/N;
dy = L/N;



F_domen_cub_2d_arr = @(A, B, arr)   A*(arr.*arr) + B*(arr.*arr.*arr.*arr)


F_domen_cub_2d = @(A, B, arr) sum(sum(F_domen_cub_2d_arr(A, B, arr)));

b_arr_new = F_domen_cub_2d_arr(A, B, b_arr)

F = F_domen_cub_2d(A, B, b_arr)

domen_optfunc = @(arr) F_domen_cub_2d(A, B, arr);

options = optimoptions(@fmincon,'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations', 3000);
% 
tic
[a_arr1, fmin]= fmincon(domen_optfunc,b_arr_0,[],[],[],[],[],[],[],options)
toc
arr_itog = reshape(a_arr1, N, N, 2)

% Очень странно работает. При начальном нуле не работает. При разных
% начальных значениях может выдавать различные доменные варианты.
% (Странности при 1, 2, 5, 0.002, 0.008



%% Interaction 3D
clc 
clear

N=7; %Кол-во атомов в стержне

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
 




% options = optimoptions(@fmincon, 'Algorithm', 'sqp','Display','iter', 'MaxFunctionEvaluations', 10^7, 'MaxIterations', 10^3, 'UseParallel', true);
% options = optimoptions(@fmincon,'Algorithm', 'sqp','Display','iter', 'MaxIterations', 1000);
options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations', 100, 'Algorithm', 'interior-point');

% 
tic
[a_arr1, fmin1]= fmincon(optfunc,work_arr,[],[],Aeq,beq, [],[],[],options);
toc


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
 




% options = optimoptions(@fmincon, 'Algorithm', 'sqp','Display','iter', 'MaxFunctionEvaluations', 10^7, 'MaxIterations', 10^3, 'UseParallel', true);
% options = optimoptions(@fmincon,'Algorithm', 'sqp','Display','iter', 'MaxIterations', 1000);
options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations',X, 'Algorithm', 'interior-point');

% 
tStart = tic;
[a_arr1, fmin1]= fmincon(optfunc,work_arr,[],[],Aeq,beq, [],[],[],options);
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
% f_uxx = func_u_xx(b_arr, dx);
% f_uyy = func_u_yy(b_arr, dy);
% f_uxy = func_u_xy(b_arr, dx, dy);
% 
% f_array = 1/2 * c11*(f_uxx.^2+f_uyy.^2)   + c12*(f_uxx.*f_uyy) + 2*c44*(f_uxy.^2);


test_arr = randi(5, 3, 3, 3, 3);

%%
N = 3;
% a = flip(test_arr, 3);
u_xy_0 = diff(test_arr(:,:,:,1), 1, 2);
u_yx_0 = diff(test_arr(:,:,:,2), 1, 1);
u_xy = func_u_3d(test_arr, 1, 2, 1, 1);

%%
function array = func_u_xx(a_arr, dx)
    array = zeros(numel(a_arr(:, 1, 1)), numel(a_arr(1, :, 1)));
    array(1:numel(a_arr(:, 1, 1))-1, :) = diff(a_arr(:, :, 1), 1, 1)/dx;
end

function array = func_u_yy(a_arr, dy)
    array = zeros(numel(a_arr(:, 1, 1)), numel(a_arr(1, :, 1)));
    array(:, 1:numel(a_arr(1, :, 1))-1) = diff(a_arr(:, :, 2), 1, 2)/dy;
end

function array = func_u_xy(a_arr, dx, dy)
    array = zeros(numel(a_arr(:, 1, 1)), numel(a_arr(1, :, 1)));
    array(:, 1:numel(a_arr(1, :, 1))-1) = array(:, 1:numel(a_arr(1, :, 1))-1) + diff(a_arr(:, :, 1), 1, 2)/dx;
    array(1:numel(a_arr(1, :, 1))-1, :) = array(1:numel(a_arr(1, :, 1))-1, :) + diff(a_arr(:, :, 2), 1, 1)/dy;
end

function array = elasticity_arr(c11,c12,c44,uxx,uyy,uxy)
    array = 1/2 * c11*(uxx.^2+uyy.^2)   + c12*(uxx.*uyy) + 2*c44*(uxy.^2);
end

function array = func_a_arr(a_arr, N)
    array = reshape(a_arr, [N, N, 2]);
end

function stop = domen_plot(arr, value, state)
    num = numel(arr)/2;
    xx = 1:1:num;
    subplot(1, 2, 1)
    plot(xx, arr(1:num, 1));
    subplot(1, 2, 2)
    plot(xx, arr(num+1:end, 1));
    drawnow;
    stop = false;
end

function array = func_arr_3d(a_arr, N)
    array = reshape(a_arr, [N, N, N, 3]);
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

function array = domen_gradient_arr_3d(G11, G14, G44, p11, p22, p33, p12, p21, p23, p32, p31, p13)
    array = 1/2*G11*(p11.^2 + p22.^2 + p33.^2) ...
        + G14*(p11.*p22 + p22.*p33 + p11.*p33) ...
        + G44*(p12.^2 + p21.^2 + p23.^2 + p32.^2 + p31.^2 + p13.^2);
end

% function value_energy = energy()
% 
% end
% 
% function value_energy = func_for_optimization()
% 
% end

%% Comments for the next time
% Рассмотреть оптимизацию функции domen_diff_3d 
% Возможно, выгоднее по времени выполнения, создать отдельные функции для
% каждого направления дифференцирования.
% Аналогично рассмотреть func_u_3d.
% Аналогичным решением может быть более эффективная общая функция
% дифференцирования.

% Убрать повторяемость разбиения первоначального 1D массива на доменный и
% упругий 3D/4D массивы. Возможно решение через создание дополнительной
% функции, вызываемой через анонимную с передачей в неё всех необходимых
% параметров.

% Посмотреть, можно ли где-то убрать повторяемые "дорогие" действия.
% Проверить возможность ускорения засчёт увеличения требуемой памяти
% (Перед этим каким-либо образом убедиться, что программа будет продолжать
% работать при больших матрицах с этими затратами по памяти)

% Изучить large-scale алгоритмы для оптимизации начальных измерений и
% добивки необходимой точности 


























