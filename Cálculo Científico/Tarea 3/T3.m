
%% 1E
x_ = [-2,-1,1,3];
fx_ = [2,18,14,82];

% DividedDifferences(x_,fx_)

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



% pol = @(x) 0.99891-0,0000344 *x-1.96668*x.^2-0.00004*x.^3+0.5516*x.^4
% pol = [0.5516 -0.00004 -1.96668 -0,0000344 0.99891]
% plot(x,pol(x))
% grid on

% polyval(pol,0.1)

res = DividedDifferences(xnod,ynod)
res = res(1,:)

%% 2B

x_ = [1,2,3,4];
fx_ = [2,2,1,11];

% SplineConstruction(x_,fx_)

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

% plot(xplot_,S,'red',xplot_,y,'blue')
% plot(xplot_,S,'red',xplot_,y,'blue',x_,fx_,'black')

%%

x_ = [1,2,3,4];
fx_ = [2,2,1,11];
xplot_ = 1:0.01:4;

Scf_= SplineConstruction(x_,fx_)

S=EvaluateSpline(x_,Scf_,xplot_);
Scf_(1,1)
S(1)
plot(xplot_,S,'red')



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

function r = fun(x)
r = 1/(25*x^2+1)
end

function res = pol(x)
% r = 0.99891-0,0000344 *x-1.96668*x^2-0.00004*x^3+0.5516*x^4 ;
res = -0.4161+1.9463*(x+1)-0.6869*(x+0.6)*(x+1)-0.8826*(x+0.6)*(x+1)*(x+0.2)+0.5516*(x-0.2)*(x+0.6)*(x+1)*(x+0.2)
end


