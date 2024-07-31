A = [1 3;-3 1;1 3];
B = [7 -1;1 6;-1 -2];
sum(sum(A.*B))
n = norm(A,1)
%%
A = [0 3 -3; -1 0 0; 0 2 2];
norm(A,2)^2
B = [-1 0 1; 0 -1 -1;1 -1 0]
cond(B)
%% 1(a)

% GAUSS
a = -1
b = 1
I = zeros(1,34);
f=@(x) abs(x)
for (n=2:1:35)
    [x,w] = mylegpts(n);
    m = length(x);
    I_t = 0;
    for (i=1:m)
        I_t = I_t+f(x(i))*w(i);
    end
    I(n-1) = I_t;
end

exct = 1
error_gauss = abs(exct - I) 

I = zeros(1,34);
for (n=2:1:35)
    [x,w] = mychebpts(n);
    m = length(x);
    I_t = 0;
    for (i=1:m)
        I_t = I_t+f(x(i))*w(i);
    end
    I(n-1) = I_t;
end

error_curt = abs(exct - I) 
semilogy(2:1:35,error_curt,'blue',2:1:35,error_gauss,'red')
legend('Clenshaw–Curtis','Gauss-Legendre')
title('Error |x|')
xlabel('n') 
ylabel('Error') 
% semilogy(2:1:35,error)
% title('Error x^{20}')
% xlabel('n') 
% ylabel('Error') 
%%

% CURT
I = zeros(1,34);
f=@(x) exp(x)
for (n=2:1:35)
    [x,w] = mychebpts(n);
    m = length(x);
    I_t = 0;
    for (i=1:m)
        I_t = I_t+f(x(i))*w(i);
    end
    I(n-1) = I_t;
end

exct = 2/21
error_curt = abs(exct - I) 
semilogy(2:1:35,error_curt,'blue',2:1:35,error_gauss,'red')
legend('Clenshaw–Curtis','Gauss-Legendre')
title('Error x^{20}')
xlabel('n') 
ylabel('Error') 
%% 1a-vi

I = zeros(1,34);
f=@(x) exp(x)
for (n=2:1:35)
    [x,w] = mylegpts(n);
    m = length(x);
    I_t = 0;
    for (i=1:m)
        I_t = I_t+f(x(i)+4)*w(i);
    end
    I(n-1) = I_t;
end

I

%% 1a ultima
I = zeros(1,34);
f=@(x) exp(-x^(-2))
for (n=2:1:35)
    [x,w] = mylegpts(n);
    m = length(x);
    I_t = 0;
    for (i=1:m)
        I_t = I_t+f(x(i)*1.5+3.5)*w(i);
    end
    I(n-1) = I_t*3/2;
end

I
%% 2
n = 10
x = -1:2/(n-1):1
%%

[xi2DQ,w2DQ]=GaussLeg2DQuad(10);
f=@(xi,eta) (xi.^10).*(eta.^2)
numInt=w2DQ*f(xi2DQ(:,1),xi2DQ(:,2))
exactInt=(1/(10+1))*(1^(10+1)-(-1)^(10+1))*(1/(2+1))*(1^(2+1)-(-1)^(2+1))

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