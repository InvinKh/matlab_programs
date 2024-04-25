%% Descend
clc
clear

% Count points in each direction
N=8;

% Create massive
work_arr = zeros(N, N, 2);
% Create corner points's deviation
deviation = 0.3;
work_arr(1, 1, 1) = deviation;
work_arr(N, 1, 1) = deviation;
work_arr(1, N, 1) = -deviation;
work_arr(N, N, 1) = -deviation;

% Calculate distance between atoms
L = 1;
dx = L/N;
dy = L/N;

% Constants for elastic-elastic
c11 = 27.5;
c12 = 17.9;
c44 = 5.43;

