%% 1(a) (ejemplo 7)

a = 2
b = 5
I = zeros(1,34);
f=@(x) exp(-x^(-2))
for (n=2:1:35)
    [x,w] = mylegpts(n);
    m = length(x);
    I_t = 0;
    for (i=1:m)
        I_t = I_t+f(0.5*(x(i)*(b-a)+a+b))*w(i);
    end
    I(n-1) = I_t*0.5*(b-a);
end

exct = 2.718513678645963
error_gauss = abs(exct - I) 

I = zeros(1,34);
for (n=2:1:35)
    [x,w] = mychebpts(n);
    m = length(x);
    I_t = 0;
    for (i=1:m)
        I_t = I_t+f(0.5*(x(i)*(b-a)+a+b))*w(i);
    end
    I(n-1) = I_t*0.5*(b-a);
end

error_curt = abs(exct - I) 
semilogy(2:1:35,error_curt,'blue',2:1:35,error_gauss,'red')
legend('Clenshaw–Curtis','Gauss-Legendre')
title('Error e^{-1/x^2}')
xlabel('n') 
ylabel('Error') 

%% 1bi

m = 20
n = 8000
[xi2DQ,w2DQ] = GaussLeg2DQuad(8);
f = @(xi,eta) (xi.^m).*(eta.^n);
numInt = w2DQ*f(xi2DQ(:,1),xi2DQ(:,2))
abs(exactInt(m,n)-numInt)

function [xi2DQ,w2DQ]=GaussLeg2DQuad(n1D)
[x,w] = mychebpts(n1D);
xi2DQ = zeros([2 n1D^2]);
w2DQ = zeros(1,n1D^2);
for i=1:1:n1D
    for j = 1:1:n1D
        c = (i-1)*n1D+j;
        xi2DQ(c,1) = x(i);
        xi2DQ(c,2) = x(j);
        w2DQ(c) = w(i)*w(j);
    end
end 
xi2DQ = xi2DQ(:, [1,2]);
end

%% FUNCIONES

function exactInt = exactInt(m,n)
exactInt = (1/(n+1))*(1^(n+1)-(-1)^(n+1))*(1/(m+1))*(1^(m+1)-(-1)^(m+1));
end

function [x,w] = mychebpts(n)
Ns=n-1;
c=zeros(n,2);
c(1:2:n,1)=(2./[1 1-(2:2:Ns).^2 ])'; c(2,2)=1;
f=real(ifft([c(1:n,:);c(Ns:-1:2,:)]));
w=[f(1,1); 2*f(2:Ns,1); f(n,1)]';
x=flip(Ns*f(1:n,2));
end

function [x,w] = mylegpts(n)
beta=.5./sqrt(1-(2*(1:n-1)).^(-2)); 
T=diag(beta,1)+diag(beta,-1);       
[V,D]=eig(T);                       
x=diag(D);                          
[x,i]=sort(x);                      
w=2*V(1,i).^2;                      
ii=1:floor(n/2);
x=x(ii);
w=w(ii);
if (mod(n,2))
    x=[x;0;-x(end:-1:1)];
    w=[w,2-sum(2*w),w(end:-1:1)];
else 
    x=[x;-x(end:-1:1)];
    w=[w,w(end:-1:1)];
end
end