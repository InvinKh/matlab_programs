%%
function [X,U,P]=temp6
N=3000;
data = zeros(idivide(int16(N),int16(100)), 2);

n=[10 10 6]; % размер решётки (n-1)
global A C Q G d % константы
d=1; % шаг решётки
% A=[-3.712 6.079 1.303 1.294 -1.950 -2.5 3.863 2.529 1.637 1.367]; % порядки!
A=[-3.712*10^(7) 6.079*10^(8) 1.303*10^(8) 1.294*10^(9) -1.950*10^(9) -2.5*10^(9) 3.863*10^(10) 2.529*10^(10) 1.637*10^(10) 1.367*10^(10)]; % порядки!
% C=[27.5 17.9 5.43]; G=[51 0 2]; Q=[14.2 -0.74 1.57];
C=[27.5 17.9 5.43]*10^(-10); G=[51 0 2]*10^(-11); Q=[14.2 -0.74 1.57]*10^(9);

[x,y,z]=ndgrid(d*(0:n(1)),d*(0:n(2)),d*(0:n(3))); X={x;y;z}; clear x y z
U=cell(3,1); P=U; E=cell(6,1);
for i=1:3, U{i}=zeros(size(X{1})); P{i}=zeros(size(X{1})); end
% задаём значения U, P при z=0
global fix
fix=false(size(X{1})); fix(:,:,1)=true;

% Сделаем сдвиг на 0.2 через meshgrid (суммарный сдвиг)
% x_fix = linspace(0, 0.2, n(1)+1);
% y_fix = linspace(0, 0.2, n(2)+1);
% [X_fix, Y_fix] = meshgrid(x_fix, y_fix);
% X_fix = reshape(X_fix, [numel(X_fix), 1]);
% Y_fix = reshape(Y_fix, [numel(Y_fix), 1]);

% U{1}(fix)=(X{1}(fix)-X_fix)*d; U{2}(fix)=(X{2}(fix)-Y_fix)*d;
% for i=1:3, P{i}(fix)=(rand(n(1:2)+1)-1/2)*d; end

% цикл градиентного спуска
F=energies(U,P); k0=0; oo=0;
% while true
for i =1:N
        oo=oo+1;
        [GU,GP,gmax]=gradflow; % градиенты
%         if gmax<1e-4*d, break, end % критерий остановки
        mu=min(1,1e-3*d/gmax); f=step_along(U,P,GU,GP,mu); f_=f;
        k=0;
        % определяем интервал одномерного поиска
        if f>F
            while f>F, f=f_; mu=mu/2; f_=step_along(U,P,GU,GP,mu); end
            mu=2*mu; k=k+1;
        else       
            while f<=f_, f_=f; mu=2*mu; f=step_along(U,P,GU,GP,mu); k=k+1; end
        end
        if k>k0,k0=k; end
        if f-2*f_+F == 0, break, end
        mu=mu/4*(1-2*(f_-F)/(f-2*f_+F));
        [F,U,P]=step_along(U,P,GU,GP,mu);
        if mod(oo,100)==0 
            data(idivide(int16(oo),int16(100)), 1) = oo;
            data(idivide(int16(oo),int16(100)), 2) = F;
            
        end
end
save(['data' point '.mat'],'data', 'X', 'U', 'P')
visual(X,U,P)
disp([oo k0])
end
%%
function [X,U,P]=my_temp6(X,U,P,N)
% N=100;
data = zeros(idivide(int16(N),int16(100)), 2);

n=[10 10 5]; % размер решётки (n-1)
global A C Q G d % константы
d=1; % шаг решётки
A=[-3.712 6.079 1.303 1.294 -1.950 -2.5 3.863 2.529 1.637 1.367]; % порядки!
C=[27.5 17.9 5.43]; G=[51 0 2]; Q=[14.2 -0.74 1.57];

% цикл градиентного спуска
F=energies(U,P); k0=0; oo=0;
% while true
for i =1:N
        oo=oo+1;
        [GU,GP,gmax]=gradflow; % градиенты
%         if gmax<1e-4*d, break, end % критерий остановки
        mu=min(1,1e-3*d/gmax); f=step_along(U,P,GU,GP,mu); f_=f;
        k=0;
        % определяем интервал одномерного поиска
        if f>F
            while f>F, f=f_; mu=mu/2; f_=step_along(U,P,GU,GP,mu); end
            mu=2*mu; k=k+1;
        else       
            while f<=f_, f_=f; mu=2*mu; f=step_along(U,P,GU,GP,mu); k=k+1; end
        end
        if k>k0,k0=k; end
        mu=mu/4*(1-2*(f_-F)/(f-2*f_+F));
        [F,U,P]=step_along(U,P,GU,GP,mu);
        if mod(oo,100)==0 
            data(idivide(int16(oo),int16(100)), 1) = oo;
            data(idivide(int16(oo),int16(100)), 2) = F;
            
        end
end
save('data.mat','data', 'X', 'U', 'P')
visual(X,U,P)
disp([oo k0])
end

%% ========================================================================
function [F,U,P]=step_along(U,P,GU,GP,mu)
% энергия после шага в направлении антиградиента: U -> (U - mu*dU)
for i=1:3, U{i}=U{i}-mu*GU{i}; P{i}=P{i}-mu*GP{i}; end
F=energies(U,P);
end

%% ========================================================================
function [F,FL,FC,FQ,FG]=energies(U,P0)
% вычисление энергии и её составляющих
global A C Q G d E P dP
% считаем производные на решётке и деформации
dU=derivatives(U); dP=derivatives(P0); P=P0;
i1=[1 5 9 6 3 2]; i2=[1 5 9 8 7 4]; % для нумерации по Фойгту
for i=1:6, E{i}=(dU{i1(i)}+dU{i2(i)})/(2-(i>3)); end
% квадраты и 4-е степени для ускорения счёта
q12=P{1}.^2; q22=P{2}.^2; q32=P{3}.^2; q14=q12.^2; q24=q22.^2; q34=q32.^2;
FL=A(1)*(q12+q22+q32)...
    +A(2)*(q14+q24+q34)...
    +A(3)*(q12.*q22+q22.*q32+q12.*q32)...
    +A(4)*(q12.*q14+q22.*q24+q32.*q34)...
    +A(5)*(q14.*(q22+q32)+q24.*(q12+q32)+q34.*(q12+q22))...
    +A(6)*q12.*q22.*q32...
    +A(7)*(q14.^2+q24.^2+q34.^2)...
    +A(8)*(q12.*q14.*(q22+q32)+q22.*q24.*(q12+q32)+q32.*q34.*(q12+q32))...
    +A(9)*(q14.*q24+q24.*q34+q14.*q34)...
    +A(10)*(q14.*q22.*q32+q24.*q12.*q32+q34.*q12.*q22);
FC=C(1)/2*(E{1}.^2+E{2}.^2+E{3}.^2)...
    +C(2)*(E{2}.*E{3}+E{1}.*E{3}+E{1}.*E{2})...
    +C(3)/2*(E{4}.^2+E{5}.^2+E{6}.^2);
FQ=-Q(1)*(E{1}.*q12+E{2}.*q22+E{3}.*q32)...
    -Q(2)*(E{1}.*(q22+q32)+E{2}.*(q12+q32)+E{3}.*(q12+q22))...
    -Q(3)*(E{6}.*P{1}.*P{2}+E{5}.*P{1}.*P{3}+E{4}.*P{2}.*P{3});
FG=G(1)/2*(dP{1,1}.^2+dP{2,2}.^2+dP{3,3}.^2)...
    +G(2)*(dP{1,1}.*dP{2,2}+dP{2,2}.*dP{3,3}+dP{1,1}.*dP{3,3})...
    +G(3)/2*(dP{1,2}.^2+dP{2,1}.^2+dP{2,3}.^2+dP{3,2}.^2+dP{3,1}.^2+dP{1,3}.^2);
FL=sum(FL(:)); FC=sum(FC(:)); FQ=sum(FQ(:)); FG=sum(FG(:));
F=FL+FC+FQ+FG; % * d^3 ?
end

%% ========================================================================
function D=derivatives(V)
% вычисляем производные векторного поля на решётке: внутренние узлы по
% формулам 2-го порядка, граничные - по формулам 1-го порядка
% компоненты поля V - массив cell(3,1)
global d
D=cell(3,3); for i=1:9, D{i}=zeros(size(V{1})); end
for i=1:3, V{i}=V{i}/d/2; end % сразу учитываем шаг решётки
for i=1:3
    D{i,1}(2:end-1,:,:)=V{i}(3:end,:,:)-V{i}(1:end-2,:,:);
    D{i,1}([1 end],:,:)=2*cat(1,V{i}(2,:,:)-V{i}(1,:,:),V{i}(end,:,:)-V{i}(end-1,:,:));
    D{i,2}(:,2:end-1,:)=V{i}(:,3:end,:)-V{i}(:,1:end-2,:);
    D{i,2}(:,[1 end],:)=2*cat(2,V{i}(:,2,:)-V{i}(:,1,:),V{i}(:,end,:)-V{i}(:,end-1,:));
    D{i,3}(:,:,2:end-1)=V{i}(:,:,3:end)-V{i}(:,:,1:end-2);
    D{i,3}(:,:,[1 end])=2*cat(3,V{i}(:,:,2)-V{i}(:,:,1),V{i}(:,:,end)-V{i}(:,:,end-1));
end
end

%% ========================================================================
function [GU,GP,gmax]=gradflow
% вычисляем градиент - производные энергии по (P1,P2,P3) и (U1,U2,U3) каждого узла
global A C Q G d E P dP fix
GP=cell(3,1); GU=cell(3,1);
% вспомогательные переменные для ускорения счёта
q2=cell(3,1); q4=q2; q6=q2; PP=q2; P2=0; P4=0; P6=0;
for i=1:3
    q2{i}=P{i}.^2; q4{i}=q2{i}.^2; q6{i}=q2{i}.^3;
    P2=P2+q2{i}; P4=P4+q4{i}; P6=P6+q6{i};
end
PP{1}=q2{2}.*q2{3}; PP{2}=q2{1}.*q2{3}; PP{3}=q2{1}.*q2{2};
EP{1}=E{6}.*P{2}+E{5}.*P{3}; EP{2}=E{6}.*P{1}+E{4}.*P{3}; EP{3}=E{5}.*P{1}+E{4}.*P{2};
% локальная часть
for i=1:3
    % производная энергии Ландау
    GP{i}=A(1)+2*A(2)*q2{i}+A(3)*(P2-q2{i})+3*A(4)*q4{i}...
        +A(5)*(2*q2{i}.*(P2-q2{i})+(P4-q4{i}))+A(6)*PP{i}+4*A(7)*q6{i}...
        +A(8)*(3*q4{i}.*(P2-q2{i})+P6-q6{i})+2*A(9)*q2{i}.*(P4-q4{i})...
        +A(10)*PP{i}.*(q2{i}+P2);
    GP{i}=2*P{i}.*GP{i};
    % вклад электрострикции
    GP{i}=GP{i}-2*P{i}.*(Q(2)*(E{1}+E{2}+E{3})-(Q(2)-Q(1))*E{i})-Q(3)*EP{i};
end
% нелокальная часть
g=G/2/d; c=C/2/d; q=Q/2/d; % если нужны конечные разности, то не делить на d!
for i=1:3
    per=[i:3 1:i-1]; % перестановка размерностей
    i1=per(1); i2=per(2); i3=per(3);
    % производные по компонентам вектора P
    % производные "вдоль", типа P1,1
    T=g(1)*dP{i1,i1}+g(2)*(dP{i2,i2}+dP{i3,i3});
    GP{i1}=GP{i1}+gradflow_local(T,i1);
    % производные "поперёк", типа P1,2 и P1,3
    T=g(3)*dP{i1,i2};
    GP{i1}=GP{i1}+gradflow_local(T,i2);
    T=g(3)*dP{i1,i3};
    GP{i1}=GP{i1}+gradflow_local(T,i3);
    % производные по компонентам вектора U (упругая энергия и электрострикция)
    T=c(1)*E{i1}+c(2)*(E{i2}+E{i3})-q(1)*q2{i1}-q(2)*(q2{i2}+q2{i3}); 
    GU{i1}=gradflow_local(T,i1);
    T=c(3)*E{i2+3}-q(3)*P{i1}.*P{i3}; % индексы согласованы с нумерацией Фойгта
    GU{i1}=GU{i1}+gradflow_local(T,i3);
    T=c(3)*E{i3+3}-q(3)*P{i1}.*P{i2};
    GU{i1}=GU{i1}+gradflow_local(T,i2);
    % GU{i1}(fix)=0; % GP{i1}(fix)=0; % для фиксированных узлов
end
% максимальная длина градиента U
gmax=sqrt(max(GU{1}(:).^2+GU{2}(:).^2+GU{3}(:).^2));
end
%% ========================================================================
function D=gradflow_local(T,dim)
% преобразование массива в разностный вдоль размерности dim 
% (вспомог. функция для gradflow)
D=zeros(size(T)); Z=false(size(T));
n=size(T); ind={1:n(1) 1:n(2) 1:n(3)};
% внутренние слои
ind{dim}=2:n(dim)-1; Z(ind{:})=true; ind0=Z; Z(:)=false;
ind{dim}=3:n(dim); Z(ind{:})=true; ind1=Z; Z(:)=false;
ind{dim}=1:n(dim)-2; Z(ind{:})=true; ind2=Z; Z(:)=false;
D(ind1)=T(ind0); D(ind2)=D(ind2)-T(ind0);
% граничные и приграничные слои
ind{dim}=[1 n(dim)]; Z(ind{:})=true; ind0=Z; Z(:)=false;
ind{dim}=[2 n(dim)]; Z(ind{:})=true; ind1=Z; Z(:)=false;
ind{dim}=[1 n(dim)-1]; Z(ind{:})=true; ind2=Z;
D(ind1)=D(ind1)+2*T(ind0); D(ind2)=D(ind2)-2*T(ind0);
end

%% ========================================================================
function visual(X,U,P)
% визуализация решётки с деформациями и поляризацией
dx=zeros(3,1); xmin=dx; X1=X; G=dx; G1=dx;
d=X{1}(2,1,1)-X{1}(1,1,1); % шаг решётки
for i=1:3, xmin(i)=min(X{i}(:)); dx(i)=max(X{i}(:))-xmin(i); end
% строим подложку
patch(xmin(1)+[-1 5 5 -1]/4*dx(1),xmin(2)+[-1 -1 5 5]/4*dx(2),[0 0 0 0],...
    'facecolor',0.9*[1 1 1])
% находим максимальные смещения и поляризации
umax=U{1}.^2+U{2}.^2+U{3}.^2; pmax=P{1}.^2+P{2}.^2+P{3}.^2;
umax=sqrt(max(umax(:))); pmax=sqrt(max(pmax(:)));
% нормируем векторные поля с учётом размера ячейки
for i=1:3
    U{i}=(d/umax)*U{i}; P{i}=(3/4*d/pmax)*P{i}; X1{i}=X{i}+U{i};
end
% строим исходную и деформированную решётки
for i=1:3
    x=X; x1=X1; per=[i:3 1:i-1];
    for j=1:3
        x{j}=permute(x{j},per); x{j}(end+1,:,:)=NaN;
        x1{j}=permute(x1{j},per); x1{j}(end+1,:,:)=NaN;
    end
    G(i)=line(x{1}(:),x{2}(:),x{3}(:),'color',[1 0 0 0.2],'linewidth',0.1);
    G1(i)=line(x1{1}(:),x1{2}(:),x1{3}(:),'color',[0 0 1 0.5],'linewidth',1);
end
axis equal, hold on
% строим векторы поляризации
for i=1:3, XP{i}=[X1{i}(:) X1{i}(:)+P{i}(:)]'; end
line(XP{1},XP{2},XP{3},'color',[0.9 0.4 0.1],'linewidth',2.5)
scatter3(X1{1}(:),X1{2}(:),X1{3}(:),'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .75 .75],'SizeData',35)
axis tight
end