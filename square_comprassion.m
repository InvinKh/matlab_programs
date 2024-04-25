%%
clc
clear

n=5; % нечетное

YX=cell(n);


for i=1:n
    for j=1:n
          YX{i,j}=zeros(1,8);
    end
end

% Про цикл сверху. массиве (n x n) точек(атомов), каждому присвоили 8 характеристик : a_x ,
% a_y, u_xx, u_xy, u_yy, sigma_xx, sigma_xy, sigma_yy (пока что всё равно 0).



YX{1,1}(1)=-0.1; % Задаём растяжение нижней границы (дели на 2, т.е. если хочешь чтобы растянулось на 2 пиши 1)



for i=1:n
    YX{1,i}(1)=YX{1,1}(1)-YX{1,1}(1)*(i-1)/(n-1)*2;
end
% Задали смещение нижней границы


for i=1:n
 YX{i,1}(6)=0;
end
for i=1:n
 YX{5,i}(8)=0;
end
for i=1:n
 YX{i,5}(6)=0;
end
for i=1:n
 YX{1,i}(8)=0;
end
%В циклах выше задаются г.у. на сигму на 4х гранях квадрата (пока что по условию везде 0, но на всякий случай оставлю)

for i=2:n-1
    for j=1:n-1
        YX{i,j}(3) = YX{i,j}(3)+(YX{i,j+1}(1)-YX{i,j}(1))/2;
        YX{i,j+1}(3) = YX{i,j+1}(3)+(YX{i,j+1}(1)-YX{i,j}(1))/2;


        YX{i,j}(5)= YX{i,j+1}(5)+(YX{i+1,j}(2)-YX{i,j}(2))/2;
        YX{i+1,j}(5)= YX{i+1,j}(5)+(YX{i+1,j}(2)-YX{i,j}(2))/2;
        YX{i,j}(5)= YX{i,j+1}(5)+(YX{i-1,j}(2)-YX{i,j}(2))/2;
        YX{i-1,j}(5)= YX{i-1,j}(5)+(YX{i-1,j}(2)-YX{i,j}(2))/2;

        YX{i,j}(4)= YX{i,j}(4)+0.5*(YX{i+1, j}(1)-YX{i, j}(1) + YX{i, j+1}(2)-YX{i, j}(2))/4;
        YX{i+1,j}(4)= YX{i+1,j}(4)+0.5*(YX{i+1, j}(1)-YX{i, j}(1) + YX{i, j+1}(2)-YX{i, j}(2))/4;
        YX{i,j}(4)= YX{i,j}(4)+0.5*(YX{i-1, j}(1)-YX{i, j}(1) + YX{i, j+1}(2)-YX{i, j}(2))/4;
        YX{i-1,j}(4)= YX{i-1,j}(4)+0.5*(YX{i-1, j}(1)-YX{i, j}(1) + YX{i, j+1}(2)-YX{i, j}(2))/4;
        YX{i,j}(4)= YX{i,j}(4)+0.5*(YX{i+1, j}(1)-YX{i, j}(1) + YX{i, j+1}(2)-YX{i, j}(2))/4;
        YX{i,j+1}(4)= YX{i,j+1}(4)+0.5*(YX{i+1, j}(1)-YX{i, j}(1) + YX{i, j+1}(2)-YX{i, j}(2))/4;
    end

end

%% задача про растяжение/сжатие стержня на определенную величину
clc 
clear

n=5; %количество атомов в столбце/строчке

YX=cell(n);


for i=1:n
    for j=1:n
          YX{i,j}=zeros(1,8);
    end
end

z=1;% Межатомное расстояние
a10= 11; %Граничные условия (на сколько оттягиваем левый конец);
a1=-11; %Гр.усл. левый конец закреплен

F= @(y)K*((y(2)-a(1))^2+(y(3)-y(2))^2+(y(4)-y(3))^2+(y(5)-y(4))^2+(y(6)-y(5))^2+(y(7)-y(6))^2+(y(8)-y(7))^2+(y(9)-y(8))^2+(a(10)-y(9))^2);

% F= @(y)(1/2*u_xx * u_xx + 1/2*u_yy *_u_yy + 2*u_yy *_u_xy + 2*u_xx *_u_xy + 2*u_xy *_u_xy + u_xx *_u_yy)

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

%%
clear
close all

% Размеры трехмерного кристалла
N = 50; % количество узлов в каждом направлении
L = 10; % длина кристалла в каждом направлении

% Координаты узлов трехмерной сетки
x = linspace(0, L, N);
y = linspace(0, L, N);
z = linspace(0, L, N);
[X, Y, Z] = meshgrid(x, y, z);

% Константы упругости
E = 1.0; % модуль Юнга
nu = 0.3; % коэффициент Пуассона

% Расчет жесткостных матриц для узлов сетки
D = E / ((1 + nu) * (1 - 2 * nu)) * ...
    [1 - nu, nu,  nu, 0,           0,  0;     nu,     1 - nu, nu, 0,           0,  0;     nu,     nu,  1 - nu, 0,           0,  0;     0,      0,  0,      (1 - 2 * nu)/2,  0,  0;     0,      0,  0,      0,           (1 - 2 * nu)/2,  0;     0,      0,  0,      0,           0,  (1 - 2 * nu)/2];

% Расчет деформации в каждом узле
strain_x = sin(X/L*2*pi); % деформация по оси X
strain_y = sin(Y/L*2*pi); % деформация по оси Y
strain_z = sin(Z/L*2*pi); % деформация по оси Z

% Расчет напряжений в каждом узле
stress_x = D(1,1) * strain_x + D(1,2) * strain_y + D(1,3) * strain_z;
stress_y = D(2,1) * strain_x + D(2,2) * strain_y + D(2,3) * strain_z;
stress_z = D(3,1) * strain_x + D(3,2) * strain_y + D(3,3) * strain_z;

% Визуализация распределения напряжений
figure;
slice(X, Y, Z, stress_x, L/2, L/2, L/2);
title('Напряжение по оси X');
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar;

figure;
slice(X, Y, Z, stress_y, L/2, L/2, L/2);
title('Напряжение по оси Y');
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar;

figure;
slice(X, Y, Z, stress_z, L/2, L/2, L/2);
title('Напряжение по оси Z');
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar;

%%
clc
clear

% Задание параметров модели
L = 10; % Размерность модели (длина ребра куба)
N = 5; % Количество атомов вдоль одной из осей (N^3 - общее количество атомов)

% Создание начальной конфигурации случайно расположенных атомов в модели
xyz = randn(N^3, 3); % Координаты атомов (случайно распределены)
xyz = xyz * L; % Масштабирование координат на размер модели

% Расчет сил действующих на каждый атом от остальных атомов
forces = zeros(N^3, 3); % Инициализация массива сил
for i = 1:N^3
    for j = 1:N^3
        if i~=j
            rij = xyz(i,:) - xyz(j,:); % Вектор между i-ым и j-ым атомами
            fij = rij ./ norm(rij)^3; % Закон Гука для силы между атомами (скачок потенциала)
            forces(i,:) = forces(i,:) - fij; % Накопление сил
        end
    end
end

% Расчет суммарного давления на атомы (вектор давлений)
pressure = zeros(N^3, 1); % Инициализация массива давлений
for i = 1:N^3
    for j = 1:N^3
        if i~=j
            rij = xyz(i,:) - xyz(j,:); % Вектор между i-ым и j-ым атомами
            fij = rij ./ norm(rij)^3; % Закон Гука для силы между атомами (скачок потенциала)
            pressure(i) = pressure(i) + dot(fij, rij); % Накопление давлений
        end
    end
end

% Вывод результатов
disp('Координаты атомов:');
disp(xyz);
disp('Силы, действующие на атомы:');
disp(forces);
disp('Давление на каждый атом:');
disp(pressure);