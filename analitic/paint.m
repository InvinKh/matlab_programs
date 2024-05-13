load("C:\Users\User\Labs\Matlab\data.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls\data_45_30_30_6.mat")

% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_10.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_20.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_30.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_40.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_45.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_50.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_60.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_70.mat")
% load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_80.mat")

%%
layer = 9;
figure
subplot(2, 2, 1)
surf(X{1}(:, :, layer), X{2}(:, :, layer), E{1}(:, :, layer))
xlabel('X') 
ylabel('Y')
title('Uxx')
subplot(2, 2, 2)
surf(X{1}(:, :, layer), X{2}(:, :, layer), E{2}(:, :, layer))
xlabel('X') 
ylabel('Y')
title('Uyy')
subplot(2, 2, 3)
surf(X{1}(:, :, layer), X{2}(:, :, layer), E{6}(:, :, layer))
xlabel('X') 
ylabel('Y')
title('Uxy')

%%
layer = 6;
figure
subplot(2, 2, 1)
surf(X{1}(:, :, layer), X{2}(:, :, layer), P{1}(:, :, layer))
xlabel('X') 
ylabel('Y')
title('Uxx')
subplot(2, 2, 2)
surf(X{1}(:, :, layer), X{2}(:, :, layer), P{2}(:, :, layer))
xlabel('X') 
ylabel('Y')
title('Uyy')
subplot(2, 2, 3)
surf(X{1}(:, :, layer), X{2}(:, :, layer), P{3}(:, :, layer))
xlabel('X') 
ylabel('Y')
title('Uxy')

%% 
% % Суммарная энергия
% figure
% surf(X{1}(:, :, 1), X{2}(:, :, 1), F(:, :, 1))
% xlabel('X') 
% ylabel('Y')
% title('F')
layer = 6;
figure
subplot(2, 2, 1)
surf(X{1}(:, :, layer), X{2}(:, :, layer), FC(:, :, layer))
xlabel('X') 
ylabel('Y')
title('FC')

subplot(2, 2, 2)
surf(X{1}(:, :, layer), X{2}(:, :, layer), FL(:, :, layer))
xlabel('X') 
ylabel('Y')
title('FL')

subplot(2, 2, 3)
surf(X{1}(:, :, layer), X{2}(:, :, layer), FQ(:, :, layer))
xlabel('X') 
ylabel('Y')
title('FQ')

subplot(2, 2, 4)
surf(X{1}(:, :, layer), X{2}(:, :, layer), FG(:, :, layer))
xlabel('X') 
ylabel('Y')
title('FG')


% Все энергии на 1 графике
figure
hold on
surf(X{1}(:, :, layer), X{2}(:, :, layer), FC(:, :, layer))
surf(X{1}(:, :, layer), X{2}(:, :, layer), FL(:, :, layer))
surf(X{1}(:, :, layer), X{2}(:, :, layer), FQ(:, :, layer))
surf(X{1}(:, :, layer), X{2}(:, :, layer), FG(:, :, layer))
legend('FC', 'FL', 'FQ', 'FG')


% figure
% surf(X{1}(:, :, 1), X{2}(:, :, 1), FC(:, :, 1))
% xlabel('X') 
% ylabel('Y')
% title('FC')
% 
% figure
% surf(X{1}(:, :, 1), X{2}(:, :, 1), FL(:, :, 1))
% xlabel('X') 
% ylabel('Y')
% title('FL')
% 
% figure
% surf(X{1}(:, :, 1), X{2}(:, :, 1), FQ(:, :, 1))
% xlabel('X') 
% ylabel('Y')
% title('FQ')
% 
% figure
% surf(X{1}(:, :, 1), X{2}(:, :, 1), FG(:, :, 1))
% xlabel('X') 
% ylabel('Y')
% title('FG')
%%
F_sum = [];
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_10.mat")
F_sum(1) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_20.mat")
F_sum(2) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_30.mat")
F_sum(3) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_40.mat")
F_sum(4) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_45.mat")
F_sum(5) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_50.mat")
F_sum(6) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_60.mat")
F_sum(7) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_70.mat")
F_sum(8) = sum(F, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_80.mat")
F_sum(9) = sum(F, "all");

%%
FG_sum = [];
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_10.mat")
FG_sum(1) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_20.mat")
FG_sum(2) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_30.mat")
FG_sum(3) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_40.mat")
FG_sum(4) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_45.mat")
FG_sum(5) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_50.mat")
FG_sum(6) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_60.mat")
FG_sum(7) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_70.mat")
FG_sum(8) = sum(FG, "all");
load("C:\Users\User\Labs\Matlab\resilience\analitic\datas\walls_without_epi\data_80.mat")
FG_sum(9) = sum(FG, "all");