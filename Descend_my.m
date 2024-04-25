%% Finish descend (adjust step in the function from 2 to 0.2 or something)
x = descend(@obj_large, x)
norm(x)

%%
clc 
clear

a10= 0.2; %Граничные условия (на сколько оттягиваем левый конец);
a1=0; %Гр.усл. левый конец закреплен

a=zeros(1,n);
a(1)=a1;
a(10)=a10;

c11 = 1;
c12 = 0.5;
c44 = 0.5;

F_cub_2d_arr = @(c11,c12,c44,uxx,uyy,uxy)   1/2 * c11*(uxx.^2+uyy.^2)   + c12*(uxx.*uyy) + 2*c44*(uxy.^2);
F_cub_2d = @(c11,c12,c44,uxx,uyy,uxy) sum(sum(F_cub_2d_arr(c11,c12,c44,uxx,uyy,uxy)));

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

function F = obj_large(x)
F = 0.1*x(1)^2 + sum( x(2:end).^2 )  ;
end



