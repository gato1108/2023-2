%% REVISIÓN

n = 10
Z = randn(n,n); % Matriz aleatoria
A = Z'*Z;  % Simétrica
b = randn(n,1);
x_real = A\b
[x,iter]=SOR(A,b,1,10^(-6))
norm(x-x_real)

%% B
A = [2 1 1 1; 4 3 3 1; 8 7 9 6; 4 7 7 8]
D = diag(diag(A));
U = triu(A);
L = tril(A);
T = eye(4)-inv(D)*A;
%%
b=rand(length(A),1);
radio_espectral = max(abs(eig(T))) % No convergeeeee
omega =  0.01:0.01:2;
re = zeros(length(omega),1);
it = zeros(length(omega),1);
for i=omega
    M = D*(1/i)+ L;
    T = eye(4)-inv(M)*A;
    re(round(i*100)) = max(abs(eig(T)));
    [x,iter] = SOR(A,b,i,10^(-8));
    it(round(i*100)) = iter;
end
%%
min(it)
min(re)
omega
%%
plot(omega,it)
xlabel('\omega')
ylabel('Iteración')
%%
plot(omega,re)
xlabel('\omega')
ylabel('Radio Espectral')
%% FUNCIÓN 

function [x,iter]=SOR(A,b,omega,tol)
D = diag(diag(A));
U = triu(A);
L = tril(A);
M = D*(1/omega)+ L;
N = M-A;
x = [zeros(length(A),1)]
Seguir = true;
iter = 0;
while Seguir
    iter = iter+1;
    x = inv(M)*(N*x+b);
    if iter > 10000
        Seguir = false;
    end
    if norm(A*x-b)<tol*norm(b)
        Seguir = false;
    end
end
end