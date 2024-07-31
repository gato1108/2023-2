%% 3A
t_=0:0.01:1
v0_=[0,0];
v1_=[1,0];
c0_=[2,1];
c1_=[1,0]+[-2,1]

cor = BezierCurve(v0_,c0_,c1_,v1_,t_)

plot(cor(1,:),cor(2,:))

%% HERMITE

% xh = [0,1];
% fh = [0,0];
% dfh = [0.5,-0.5];

% Her = Hermite(xh,fh,dfh,t_);

% plot(t_,Her,'blue',cor(1,:),cor(2,:),'red')

xh = [0,1];
fh = [0,1];
dfh = [6,6];

Her_x = Hermite(xh,fh,dfh,t_);

xh = [0,1];
fh = [0,0];
dfh = [3,-3];

Her_y = Hermite(xh,fh,dfh,t_);

plot(Her_x,Her_y,'blue',cor(1,:),cor(2,:),'red')

%%
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

%% FUNCIONES
function v=BezierCurve(v0,c0,c1,v1,t)
n = length(t)-1
v = zeros(2,n+1)
for i=1:(n+1)
    tt = t(i);
    v(:,i) = (1-tt)^3*v0+3*(1-tt)^2*tt*c0+3*(1-tt)*tt^2*c1+tt^3*v1;
end
end

function H=Hermite(xk,f,df,x)
% vector xk: [x_0,...,x_n]
% vector f: [f(x_0),...,f(x_n)]
% vector df: [f'(x_0),...,f'(x_n)]
% x: discrete data points (to evaluate) 
n=length(xk)-1; % number of interpolating points (minus 1)
m=length(x);    % number of discrete data points
Lk=ones(n+1,m); % Lagrange basis polynomials
a=zeros(n+1,m); % basis polynomials alpha(x)
b=zeros(n+1,m); % basis polynomials beta(x)    
H=zeros(1,m);   % Hermite interpolation polynomial Q(x)
for k=1:n+1 %from 0 to n (shifted by 1 for indexing)
   dL=0;        % derivative of k-th Lagrange basis at x_k
   for j=1:n+1 %from 0 to n (shifted by 1 for indexing)
       if j~=k
           dL=dL+1/(xk(k)-xk(j));
           Lk(k,:)=Lk(k,:).*(x-xk(j))/(xk(k)-xk(j)); 
       end
   end
   Lsq=Lk(k,:).^2;
   b(k,:)=(x-xk(k)).*Lsq;           % basis polynomial alpha(x)
   a(k,:)=(1-2*(x-xk(k))*dL).*Lsq;  % basis polynomial beta(x)
   H=H+a(k,:)*f(k)+b(k,:)*df(k);    % Hermite polynomial H(x)
end
end