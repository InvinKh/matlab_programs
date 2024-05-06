scale = 10^(-10);
global A C Q G d % константы
C=[27.5 17.9 5.43]*10^(9); Q=[14.2 -0.74 1.57]*10^(7);
P=[0 0.1616];
Uxx=linspace(-.0004, .0004, 50);
Uyy=linspace(-.0004, .0004, 50)';
% Uyy=linspace(0, 0, 1000);
[X, Y] = meshgrid(Uxx, Uyy);
Uxy=0;
F=zeros(numel(Uxx), numel(Uyy));
for i=1:numel(Uxx)
    for j=1:numel(Uyy)
        F(i, j) = energies([Uxx(i) Uyy(j) Uxy], P);
    end
end
% surf(X, Y, F)
contour(X, Y, F, 100)
xlabel('Uxx') 
ylabel('Uyy') 
grid on
a = min(min(F));
ind = find(F==a);



function [F,FQ]=energies(U,P0)
% вычисление энергии и её составляющих
global C Q d E P dP
% квадраты и 4-е степени для ускорения счёта
q12=P0(1)^2; q22=P0(2)^2;

FC=C(1)/2*(U(1)^2+U(2)^2)...
    +C(2)*(U(1)*U(2))...
    +C(3)/2*(U(3)^2);
FQ=-Q(1)*(U(1)*q12+U(2).*q22)...
    -Q(2)*(U(1).*q22+U(2).*q12)...
    -Q(3)*(U(3).*P0(1).*P0(2));
% if Q(1)*(U(1)*q12+U(2).*q22)~=0
%     disp(1)
% end
% if Q(3)*(U(3).*P0(1).*P0(1))~=0
%     disp(2) 
% end


% FQ=sum(FQ(:))
% F=FQ;
% F=FC;
F=FC+FQ; % * d^3 ?
% F=F/numel(P0(1));
end