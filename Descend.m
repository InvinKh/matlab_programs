%% TRY FMINCON
options = optimoptions('fmincon','display','iter','MaxFunctionEvaluations',1e6)
x = fmincon(@obj_large, ones(1,50000), [],[],[],[],[],[],[], options);
norm(x)

%% TRY DESCEND WITH THE SAME FUNCTION (use step_mag=2 or larger).
x = descend(@obj_large, ones(1,50000))
%x = descend(@obj_large, x)
norm(x)

%% Finish descend (adjust step in the function from 2 to 0.2 or something)
x = descend(@obj_large, x)
norm(x)

%%

function x_cur = descend(F,x0)
x_cur = x0;
dxmag = 0.00001;
step_mag = 2; % Adjust as you like
for i_descend=1:20 % Adjust as you like
    g = grad(F,x_cur,dxmag);
    g_dir = g/norm(g);
    dx = -g_dir*step_mag;
    x_cur = x_cur+dx;
    
    F_cur = F(x_cur);
    %disp([i_descend, F_cur, x_cur]);
    disp([i_descend, F_cur])
end
end


function g=grad(F,x0,dxmag)
g = 0*x0;
for i=1:numel(x0)
    F0 = F(x0);
    x1 = x0;
    x1(i) = x1(i) + dxmag;
    F1 = F(x1);
    g(i) = (F1-F0)/dxmag;
end
end

function F = obj_small(x)
F = 0.1*x(1)^2 + x(2)^2;
end

function F = obj_large(x)
F = 0.1*x(1)^2 + sum( x(2:end).^2 )  ;
end



