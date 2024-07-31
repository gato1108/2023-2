%% 1E

xnod = -1:(2/5):1
ynod = [0,0,0,0,0,0];
for i = 1:6
    ynod(i) = cos(2*xnod(i)) ;
end
x = -1:.01:1;
y_ = zeros(1,length(x));
for i = 1:(length(x))
    y_(i) = pol(x(i));
end

plot(x,y_,x,cos(2*x))
grid on
legend('Polinomio','cos(2x)')

%% 2B

x_ = [-1,-0.6,-0.2,0.2,0.6,1];
fx_ = [0,0,0,0,0,0];
for i=1:6
    fx_(i) = 1/(1+25*x_(i)^2);
end
Scf_ = SplineConstruction(x_,fx_);
xplot_ = -1:0.01:1;
y = zeros(1,201)
for i=1:201
    y(i)= fun(xplot_(i));
end


S = EvaluateSpline(x_,Scf_,xplot_);

plot(xplot_,S,'blue',xplot_,y,'red')
legend('Aproximación Spline','Función Original')

%% 3A

t_=0:0.01:1
v0_=[0,0];
v1_=[1,0];
c0_=[2,1];
c1_=[1,0]+[-2,1]

cor = BezierCurve(v0_,c0_,c1_,v1_,t_)

xh = [0,1];
fh = [0,1];
dfh = [6,6];

Her_x = Hermite(xh,fh,dfh,t_);

xh = [0,1];
fh = [0,0];
dfh = [3,-3];

Her_y = Hermite(xh,fh,dfh,t_);

plot(cor(1,:),cor(2,:))
subplot(2,2,1);
plot(cor(1,:),cor(2,:),'red')
title('Bezier')

subplot(2,2,2);
plot(Her_x,Her_y,'blue')
title('Hermite')

%% 3B

R = 1
H = 4/3*(sqrt(2)-1)
t_=0:0.01:1;
v0_=[R;0];
v1_=[0;R];
c0_=[R;H];
c1_=[H;R]

cor = BezierCurve(v0_,c0_,c1_,v1_,t_);

t=linspace(0,1);
circ_x = cos(t*pi/2);
circ_y = sin(t*pi/2);

plot(cor(1,:),cor(2,:),'red',circ_x,circ_y,'blue')
legend('Bezier','Círculo Real')

mas_lejano = (cor(:,1)+cor(:,2))/2
error = abs(1-norm(mas_lejano))

%% 3C
load('C:\Users\Panorámica\Desktop\2023 S2\Cálculo Científico\Tarea 3\AttachedFiles\Glyph-d.mat')
t=0:0.01:1;
figure
hold on
for j=1:size(V0,2);
    v=BezierCurve(V0(:,j),C0(:,j),C1(:,j),V1(:,j),t);
    plot(v(1,:),v(2,:),'b')
    axis equal
end

load('C:\Users\Panorámica\Desktop\2023 S2\Cálculo Científico\Tarea 3\AttachedFiles\Glyph-p.mat')
for j=1:size(V0,2);
    v=BezierCurve(V0(:,j),C0(:,j),C1(:,j),V1(:,j),t);
    v(1,:) = (v(1,:)+ones(1,101)*500);
    plot(v(1,:),v(2,:),'b')
    axis equal 
end

%% FUNCIONES 

function tab = DividedDifferences(x,fx)
n = length(x);
tab = zeros(n,n+1);
tab(:,1) = x ;
tab(:,2) = fx;
for j = 3:n+1
    for i = 1:(n+2-j)
        tab(i,j) = (tab(1+i,j-1)-tab(i,j-1))/(tab(i+j-2,1)-tab(i,1));
    end
end
end

function res = pol(x)
res = -0.4161+1.9463*(x+1)-0.6869*(x+0.6)*(x+1)-0.8826*(x+0.6)*(x+1)*(x+0.2)+0.5516*(x-0.2)*(x+0.6)*(x+1)*(x+0.2)
end

function S=EvaluateSpline(x,Scf,xplot)
n = length(x)-1;
t_len = length(xplot);
S = ones(1,t_len);
for i = 1:t_len
    t = xplot(i);
    for j = 1:n
        if x(j) <= t & x(j+1) >= t
            d = t-x(j);
            S(i) = Scf(j,1)+Scf(j,2)*d+Scf(j,3)*d^2+Scf(j,4)*d^3;
        end
    end
end
end


function res=SplineConstruction(x,fx)
n = length(x)-1;
h = zeros(1,n);
for i = 1:(n)
    h(i) = x(i+1)-x(i);
end
mat = zeros(n+1,n+1);
b_ = zeros(1,n+1);
for i=1:(n+1)
    if i==1
        mat(i,i) = 1;
    elseif i==(n+1)
        mat(i,i) = 1;
    else 
        mat(i,i) = 2*(h(i)+h(i-1));
        mat(i,i-1) = h(i-1);
        mat(i,1+i) = h(i);
        b_(i) = 3*((fx(i+1)-fx(i))/h(i)-(fx(i)-fx(i-1))/h(i-1));
    end
end
c = linsolve(mat,b_.');
b = zeros(n,1);
d = zeros(n,1);

for i=1:n
    b(i) = (fx(i+1)-fx(i))/h(i)-h(i)*(2*c(i)+c(i+1))/3;
    d(i) = (c(i+1)-c(i))/(3*h(i))

res = zeros(n,4) 
res(:,1) = fx(1:n)
res(:,2) = b
res(:,3) = c(1:n)
res(:,4) = d
end
end

function r = fun(x)
r = 1/(25*x^2+1)
end

function v=BezierCurve(v0,c0,c1,v1,t)
n = length(t)-1
v = zeros(2,n+1)
for i=1:(n+1)
    tt = t(i);
    v(:,i) = (1-tt)^3*v0+3*(1-tt)^2*tt*c0+3*(1-tt)*tt^2*c1+tt^3*v1;
end
end

function H=Hermite(xk,f,df,x)  
n=length(xk)-1; 
m=length(x);   
Lk=ones(n+1,m);
a=zeros(n+1,m);
b=zeros(n+1,m); 
H=zeros(1,m);  
for k=1:n+1 
   dL=0;       
   for j=1:n+1 
       if j~=k
           dL=dL+1/(xk(k)-xk(j));
           Lk(k,:)=Lk(k,:).*(x-xk(j))/(xk(k)-xk(j)); 
       end
   end
   Lsq=Lk(k,:).^2;
   b(k,:)=(x-xk(k)).*Lsq;          
   a(k,:)=(1-2*(x-xk(k))*dL).*Lsq;  
   H=H+a(k,:)*f(k)+b(k,:)*df(k);    
end
end