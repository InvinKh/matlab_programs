%%
n=[100 12 2]; % размер решётки (n-1)
scale = 10^(-10);
d=10^(-9); % шаг решётки
d = d/scale;
[x,y,z]=ndgrid(d*(0:n(1)),d*(0:n(2)),d*(0:n(3))); X={x;y;z}; clear x y z
dri=700;
part = false(size(X{1}));
ang=45;
ang=ang/180*pi;
part(find(tan(ang)*(X{1})>=X{2}+dri))=true;

%%
N=2000;
% data = zeros(idivide(int16(N),int16(100)), 2);
data = [];

% перейдём в ангстремы
global scale
scale = 10^(-10);

n=[100 12 2]; % размер решётки (n-1)
global A C Q G d % константы
d=10^(-9); % шаг решётки
d = d/scale;
A=[-3.712*10^(7) 6.079*10^(8) 1.303*10^(8)*40 1.294*10^(9) -1.950*10^(9) -2.5*10^(9) 3.863*10^(10) 2.529*10^(-10) 1.637*10^(-10) 1.367*10^(-10)]; % порядки!
% A=[-3.712*10^(7) 6.079*10^(8) 1.303*10^(8) 1.294*10^(9) -1.950*10^(9) -2.5*10^(9) 3.863*10^(10) 2.529*10^(-10) 1.637*10^(-10) 1.367*10^(-10)]; % порядки!

C=[27.5 17.9 5.43]*10^(9); 
G=[51 0 100]*10^(-11)/scale^2*1;
Q=[-6.236 -4.059 1.57]*10^(9)*0.9;
% Q=[0.795 -1.222 1.57]*10^(9);
% Q=[14.2 -0.74 1.57]*10^(9);
global E
[x,y,z]=ndgrid(d*(0:n(1)),d*(0:n(2)),d*(0:n(3))); X={x;y;z}; clear x y z
U=cell(3,1); P=U; E=cell(6,1);
for i=1:3, U{i}=zeros(size(X{1})); P{i}=zeros(size(X{1})); end
% задаём значения U, P при z=0
global fix
fix=false(size(X{1})); fix(:,:,1)=true;

del = -5;
del = del/100;
x_fix = linspace(-del*d/2, del*d/2, n(1)+1);
y_fix = linspace(-del*d/2, del*d/2, n(2)+1);
[X_fix, Y_fix] = meshgrid(x_fix, y_fix);

U{1}(fix)=(X_fix)';
U{2}(fix)=(Y_fix)';
% U{1}(fix)=0;
% U{2}(fix)=0;
U{3}(fix)=0;
visual(X,U,P)

%% ========================================================================
function visual(X,U,P)
figure
% визуализация решётки с деформациями и поляризацией
dx=zeros(3,1); xmin=dx; X1=X; G=dx; G1=dx;
d=X{1}(2,1,1)-X{1}(1,1,1); % шаг решётки
for i=1:3, xmin(i)=min(X{i}(:)); dx(i)=max(X{i}(:))-xmin(i); end
% % строим подложку
% patch(xmin(1)+[-1 5 5 -1]/4*dx(1),xmin(2)+[-1 -1 5 5]/4*dx(2),[0 0 0 0],...
%     'facecolor',0.9*[1 1 1])
% находим максимальные смещения и поляризации
umax=U{1}.^2+U{2}.^2+U{3}.^2; pmax=P{1}.^2+P{2}.^2+P{3}.^2;
umax=sqrt(max(umax(:))); pmax=sqrt(max(pmax(:)));
% нормируем векторные поля с учётом размера ячейки
for i=1:3
    X1{i}=X{i}+U{i};
    U{i}=(d/umax)*U{i}; 
    P{i}=(3/4*d/pmax)*P{i}; 
%     X1{i}=X{i}+U{i};
end

% строим исходную и деформированную решётки
subplot(2, 1, 1)
axis equal, hold on
% for i=1:3
%     x=X; x1=X1; per=[i:3 1:i-1];
%     for j=1:3
%         x{j}=permute(x{j},per); x{j}(end+1,:,:)=NaN;
%         x1{j}=permute(x1{j},per); x1{j}(end+1,:,:)=NaN;
%     end
%     G(i)=line(x{1}(:),x{2}(:),x{3}(:),'color',[1 0 0 0.2],'linewidth',0.1);
%     G1(i)=line(x1{1}(:),x1{2}(:),x1{3}(:),'color',[0 0 1 0.5],'linewidth',1);
% end

% for i=1:3, XU{i}=[X{i}(:) X{i}(:)+U{i}(:)]'; end
for i=1:3, XU{i}=[X{i}(:) X{i}(:)+U{i}(:)]'; end
line(XU{1},XU{2},XU{3},'color',[0.9 0.4 0.1],'linewidth',2.5)
scatter3(X{1}(:),X{2}(:),X{3}(:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75],'SizeData', 15)

% строим векторы поляризации
subplot(2, 1, 2)
axis equal, hold on
for i=1:3, XP{i}=[X1{i}(:) X1{i}(:)+P{i}(:)]'; end
line(XP{1},XP{2},XP{3},'color',[0.9 0.4 0.1],'linewidth',2.5)
scatter3(X1{1}(:),X1{2}(:),X1{3}(:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75],'SizeData',15)
axis tight
end

