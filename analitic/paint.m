load("C:\Users\User\Labs\Matlab\data.mat")
figure
surf(X{1}(:, :, 1), X{2}(:, :, 1), F(:, :, 1))
xlabel('X') 
ylabel('Y')
title('F')

figure
subplot(2, 2, 1)
surf(X{1}(:, :, 1), X{2}(:, :, 1), FC(:, :, 1))
xlabel('X') 
ylabel('Y')
title('FC')

subplot(2, 2, 2)
surf(X{1}(:, :, 1), X{2}(:, :, 1), FL(:, :, 1))
xlabel('X') 
ylabel('Y')
title('FL')

subplot(2, 2, 3)
surf(X{1}(:, :, 1), X{2}(:, :, 1), FQ(:, :, 1))
xlabel('X') 
ylabel('Y')
title('FQ')

subplot(2, 2, 4)
surf(X{1}(:, :, 1), X{2}(:, :, 1), FG(:, :, 1))
xlabel('X') 
ylabel('Y')
title('FG')

figure
hold on
surf(X{1}(:, :, 1), X{2}(:, :, 1), FC(:, :, 1))
surf(X{1}(:, :, 1), X{2}(:, :, 1), FL(:, :, 1))
surf(X{1}(:, :, 1), X{2}(:, :, 1), FQ(:, :, 1))
surf(X{1}(:, :, 1), X{2}(:, :, 1), FG(:, :, 1))
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
