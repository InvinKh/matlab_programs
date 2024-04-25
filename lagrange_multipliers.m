%%
clc
clear
a = 1;
b = 7;
c = 30;
my_fun = @(x, y, a, b, c) a*(x^2+y^2) + b*(x+y) + c;
my_opt_fun = @(x) my_fun(x(1), x(2), a, b, c);
opt_fun = @(x) func(x, a, b, c);
x = 0;
y = 0;
lambda = 0;
args(1) = x;
args(2) = y;
args(3) = lambda;
opt = @(args) lagr_func(args);

% opt = @(args) lagr_grad_func(args);

% options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Display','iter', 'MaxIterations', 100, 'Algorithm', 'trust-region-reflective');
% options = optimoptions(@fmincon, 'Display','iter', 'Algorithm', 'trust-region-reflective', 'SpecifyObjectiveGradient', true, 'CheckGradients', true);
% options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Display','iter', 'Algorithm', 'sqp');
% options = optimoptions(@fmincon, 'Algorithm', 'trust-region-reflective', 'SpecifyObjectiveGradient', true, 'CheckGradients', true, 'MaxIterations', 100);
options = optimoptions(@fmincon, 'MaxFunctionEvaluations',20000000, 'Algorithm', 'sqp', 'MaxIterations', 100);

tic
[a_arr1,fval,exitflag,output,lambda,grad,hessian] = fmincon(opt,args,[],[],[],[],[],[],[],options);
toc
%%
lagr_func(args);
%%
function [F, grad] = func(args, a, b, c)
    x = args(1);
    y = args(2);
    fun = @(x, y, a, b, c) a*(x^2+y^2) + b*(x+y) + c;
    F = fun(x, y, a, b, c);
    grad(1) = 2*a*x+b;
    grad(2) = 2*a*y+b;
end

function F = lagr_func(args)
    x = args(1);
    y = args(2);
    lambda = args(3);
    f = @(x, y) x^2+y^2;
    g = @(x, y) x^2+y^2-1;
    F = f(x, y) + lambda*g(x, y);
    args
end

function [F, grad] = lagr_grad_func(args)
    x = args(1);
    y = args(2);
    lambda = args(3);
    f = @(x, y) x^2+y^2;
    g = @(x, y) x^2+y^2-1;
    F = f(x, y) + lambda*g(x, y);
    grad(1) = 2*x+2*x*lambda;
    grad(2) = 2*y+2*y*lambda;
    grad(3) = g(x, y);
    args
end
