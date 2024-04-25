%%
clc
clear

x = 0;
y = 0;
lambda = 0;
% args(1) = x;
% args(2) = y;
% args(3) = lambda;

args(1) = x;
args(2) = lambda;


% opt = @(args) lagr_func(args);
% opt = @(args) lagr_grad_func(args);
% opt = @(args) lagr_grad_1D_func(args);
opt = @(args) lagr_1D_func(args);


% options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Algorithm', 'trust-region-reflective', 'SpecifyObjectiveGradient', true, 'CheckGradients', true, 'MaxIterations', 100);
options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Algorithm', 'sqp', 'MaxIterations', 100);

tic
[a_arr1,fval,exitflag,output,lambda,grad,hessian] = fmincon(opt,args,[],[],[],[],[],[],[],options);
toc
%%
[x,lambda]=meshgrid(-3:.1:2, -3:.1:2)
F = x*0;
for i = 1:numel(x)
    F(i) = lagr_1D_func([x(i) lambda(i)]);
end

figure;
contour(x, lambda, F, 50)
%%
function F = lagr_1D_func(args)
    x = args(1);
    lambda = args(2);
    f = @(x) x^2;
    g = @(x) x-1;
    F = f(x) + lambda*g(x);
end

function [F, grad] = lagr_grad_1D_func(args)
    x = args(1);
    lambda = args(2);
    f = @(x) x^2;
    g = @(x) x-1;
    F = f(x) + lambda*g(x);
    grad(1) = 2*x+lambda;
    grad(2) = g(x);
    args
end

function [F, grad] = lagr_grad_func(args)
    x = args(1);
    y = args(2);
    lambda = args(3);
    f = @(x, y) x^2+y^2;
    g = @(x, y) x-1;
    F = f(x, y) + lambda*g(x, y);
    grad(1) = 2*x+lambda;
    grad(2) = 2*y;
    grad(3) = g(x, y);
    args
end