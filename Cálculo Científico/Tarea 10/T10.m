%% 1C
tol = 10^(-6)
A = [-11 10 14;6 -2 -6;-12 10 15]
B = [14 -2 -13;-8 2 8;16 -2 -15]
eigQR(A)
eig(A)
eigQR(B)
eig(B)
%% 1D
x = [2,2,1];
v = householder(x)
(eye(3)-2*v'*v)*x'

%% 1D
A = magic(6)
prueba = hessenberg(A)
hess(A)

%% 2A
fun = @(x) 2*x.^3-4*x+4;
der = @(x) 6*x.^2-4;

x = i+1i
% x = 0.1
% x = 0.5
for i = 1:100
x = x-fun(x)/der(x);
end
x

%% 2B

ejex = -1:0.01:1
ejey = -1:0.01:1
for i = 1:201
    x = 0;
    g = @(x) x-ejex(i)*(2*x.^3-4*x+4);
    for j = 1:100
    x = g(x);
    end
    ejey(i) = x;
end

plot(ejex,ejey)

function H = hessenberg(A)
H = A
m = size(A,1)
for k = 2:m-1
    x = H(k:m,k-1);
    v = householder(x);
    Q = eye(m);
    F = (eye(m+1-k)-2*v*v');
    Q(k:m,k:m) = F;
    H = Q'*H*Q;
end
end

function v = householder(x)
n = size(x,2);
v = x;
v(1) = x(1)+sgn(x(1))*norm(x);
v = v/(norm(v));
end

function s = sgn(x)
if x<0
    s = -1;
else 
    s = 1;
end
end

function lambda = eigQR(A)
n = size(A,1)
Seguir = true;
iter = 0;
Q = eye(n);
Tn = A;
while Seguir
    Tp = Tn;
    iter = iter+1;
    V = A*Q;
    [Q,R] = qr(V);
    Tn = transpose(Q)*A*Q;
    if iter > 10000
        Seguir = false;
    end
    if norm(diag(Tn)-diag(Tp))<10^(-6)
        Seguir = false;
    end
end
lambda = diag(Tn);
end