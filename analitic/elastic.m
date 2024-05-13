function [X,U,P]=elastic
N=5000;
% data = zeros(idivide(int16(N),int16(100)), 2);
data = [];
global scale
scale = 10^(-10)

n=[20 20 5]; % размер решётки (n-1)
global C d % константы
d=10^(-8); % шаг решётки
d=d/scale
C=[27.5 17.9 5.43]*10^(9);
% d=1; % шаг решётки
% C=[27.5 17.9 5.43];

[x,y,z]=ndgrid(d*(0:n(1)),d*(0:n(2)),d*(0:n(3))); X={x;y;z}; clear x y z
U=cell(3,1);
for i=1:3, U{i}=zeros(size(X{1}));end
% задаём значения U при z=0
global fix
fix=false(size(X{1})); fix(:,:,1)=true;

% Сделаем сдвиг на 5% через meshgrid (суммарный сдвиг в процентах)
del = 50;
del = del/100;
x_fix = linspace(-del*d/2, del*d/2, n(1)+1);
y_fix = linspace(-del*d/2, del*d/2, n(2)+1);
[X_fix, Y_fix] = meshgrid(x_fix, y_fix);
% X_fix = reshape(X_fix, [numel(X_fix), 1]);
% Y_fix = reshape(Y_fix, [numel(Y_fix), 1]);

for i=1:3
    U{i}=(rand(n+1)-1/2)*d*10^(-2);
end

U{2}(fix)=(X_fix);
U{1}(fix)=(Y_fix);
U{3}(fix)=0;
% U{1}
% X{1}(:, :, 1)+U{1}(:, :, 1)

F=energies(U); 
k0=0; 
oo=0;

% while true
for i =1:N
        oo=oo+1;
        [GU,gmax]=gradflow; % градиенты
        if gmax<1e-4*d, break, end % критерий остановки
        mu=min(1,1e-3*d/gmax); f=step_along(U,GU,mu); f_=f;
        k=0;
        % определяем интервал одномерного поиска
        if f>F
            while f>F, f=f_; mu=mu/2; f_=step_along(U,GU,mu); end
            mu=2*mu; k=k+1;
        else       
            while f<=f_, f_=f; mu=2*mu; f=step_along(U,GU,mu); k=k+1; end
        end
        if k>k0,k0=k; end
        if f-2*f_+F == 0
            disp([f f_ F])
            break
        end
        mu=mu/4*(1-2*(f_-F)/(f-2*f_+F));
        [F,U]=step_along(U,GU,mu);
        if mod(oo,100)==0 
            data(idivide(int16(oo),int16(100)), 1) = oo;
            data(idivide(int16(oo),int16(100)), 2) = F;
            
        end
%         F
%         GU{1}
end
% save('data_elastic.mat','data', 'X', 'U')
% U{1}
[F, FC, dU] = arr_energies(U);
save('data_el_en.mat', "F", "dU", "U")
visual(X,U)
disp([oo k0])
end

%% ========================================================================
function [F,FC, dU]=arr_energies(U)
% function [F,FL]=energies(U,P0)
% вычисление энергии и её составляющих
global C d E
% считаем производные на решётке и деформации
dU=derivatives(U);
i1=[1 5 9 6 3 2]; i2=[1 5 9 8 7 4]; % для нумерации по Фойгту
for i=1:6, E{i}=(dU{i1(i)}+dU{i2(i)})/(2-(i>3))/d; end
% квадраты и 4-е степени для ускорения счёта
FC=C(1)/2*(E{1}.^2+E{2}.^2+E{3}.^2)...
    +C(2)*(E{2}.*E{3}+E{1}.*E{3}+E{1}.*E{2})...
    +C(3)/2*(E{4}.^2+E{5}.^2+E{6}.^2);
% FC=sum(FC(:));
F=FC; % * d^3 ?
F=F/numel(U{1});
end

%% ========================================================================
function [F,FC]=energies(U)
% function [F,FL]=energies(U,P0)
% вычисление энергии и её составляющих
global C d E
% считаем производные на решётке и деформации
dU=derivatives(U);
i1=[1 5 9 6 3 2]; i2=[1 5 9 8 7 4]; % для нумерации по Фойгту
for i=1:6, E{i}=(dU{i1(i)}+dU{i2(i)})/(2-(i>3)); end
% квадраты и 4-е степени для ускорения счёта
FC=C(1)/2*(E{1}.^2+E{2}.^2+E{3}.^2)...
    +C(2)*(E{2}.*E{3}+E{1}.*E{3}+E{1}.*E{2})...
    +C(3)/2*(E{4}.^2+E{5}.^2+E{6}.^2);
FC=sum(FC(:));
F=FC; % * d^3 ?
F=F/numel(U{1});
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
function [GU,gmax]=gradflow
% вычисляем градиент - производные энергии по (U1,U2,U3) каждого узла
global C d E fix
GU=cell(3,1);
% вспомогательные переменные для ускорения счёта
% нелокальная часть
c=C/2/d; % если нужны конечные разности, то не делить на d!
for i=1:3
    per=[i:3 1:i-1]; % перестановка размерностей
    i1=per(1); i2=per(2); i3=per(3);
    % производные по компонентам вектора U (упругая энергия и электрострикция)
    T=c(1)*E{i1}+c(2)*(E{i2}+E{i3}); 
    GU{i1}=gradflow_local(T,i1);
    T=c(3)*E{i2+3}; % индексы согласованы с нумерацией Фойгта
    GU{i1}=GU{i1}+gradflow_local(T,i3);
    T=c(3)*E{i3+3};
    GU{i1}=GU{i1}+gradflow_local(T,i2);
    GU{i1}(fix)=0;  % для фиксированных узлов
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
function visual(X,U)
% визуализация решётки с деформациями и поляризацией
dx=zeros(3,1); xmin=dx; X1=X; G=dx; G1=dx;
d=X{1}(2,1,1)-X{1}(1,1,1); % шаг решётки
for i=1:3, xmin(i)=min(X{i}(:)); dx(i)=max(X{i}(:))-xmin(i); end
% строим подложку
patch(xmin(1)+[-1 5 5 -1]/4*dx(1),xmin(2)+[-1 -1 5 5]/4*dx(2),[0 0 0 0],...
    'facecolor',0.9*[1 1 1])
% % находим максимальные смещения и поляризации
% umax=U{1}.^2+U{2}.^2+U{3}.^2;
% umax=sqrt(max(umax(:)));
% нормируем векторные поля с учётом размера ячейки
for i=1:3
    X1{i}=X{i}+U{i};
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
    xlabel('-2\pi < x < 2\pi') 
    ylabel('Sine and Cosine Values') 
end
xlabel('X, A.U.') 
ylabel('Y, A.U.') 
zlabel('Z, A.U.')
end
%% ========================================================================
function [F,U]=step_along(U,GU,mu)
% энергия после шага в направлении антиградиента: U -> (U - mu*dU)
for i=1:3, U{i}=U{i}-mu*GU{i}; end
F=energies(U);
end