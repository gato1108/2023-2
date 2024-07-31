n=200; 
Z=randn(n,n); % Matriz aleatoria
A=Z'*Z;  % Simétrica
b=randn(n,1); % Vector aleatorio
A2=A; 
A2(n,1)=A2(n,1)/2; % Dividir último elemento de primera columna por 2
I=eye(n); % Identidad
emin=min(eig(A)); % Min VP
A3=A-0.9*emin*I;
A4=A-1.1*emin*I;
A5=triu(A); % Triangular superior
A6=A5; 
A6(n,1)=A5(1,n); % Elemento en esquina
tic; for i=1:1000; x=A\b; end; toc; % Multiplicar tiempo por 1000
tic; for i=1:1000; x=A2\b; end; toc;
tic; for i=1:1000; x=A3\b; end; toc;
tic; for i=1:1000; x=A4\b; end; toc;
tic; for i=1:1000; x=A5\b; end; toc;
tic; for i=1:1000; x=A6\b; end; toc;

%% 2
Z=randn(10); A=Z'*Z; [L,U]=directlu(A); norm(A-L*U)/norm(A)

%% 3
A=randn(10); [L,U,P]=lupp(A); norm(P*A-L*U)/norm(A)

%% 4
A = [1 0 0; 2 10^(-20) 3; 0 1 2];
b = [1;-1;0];
[Ld,Ud]=directlu(A);
[Lp,Up,P]=lupp(A);
x_exc = A\b
xd=Ud\(Ld\b)
xp=Up\(Lp\P*b)
cond(Ld)
cond(Ud)
cond(Lp)
cond(Up)
cond(A)

%% 5

T = zeros(1,99);
for n = 40:20:2000
A = randn(n,n);
b = randn(n,1);
tic; x=A\b; 
T((n-20)/20)=toc;
end

log_n = log10(40:20:2000)
log_T = log10(T)
r = reg_lin(log_n(34:99),log_T(34:99));
b = r(1)
a = r(2)
x = 40:0.1:2000;
y = 10^(b)*x.^(a);
loglog(40:20:2000,T,x,y)
legend('curva', 'mejor aproximación lineal')
xlabel('log(N)')
ylabel('log(tiempo)')

%% 6
Z=randn(10); 
A=Z'*Z; 
[L,D]=myldl(A);
norm(A-L*D*L')/norm(A)

function [L,D]=myldl(A)
n = size(A); n = n(1);
L = eye(n); U = A; 
for k = 1:1:n-1
    for j = k+1:1:n
        L(j,k) = U(j,k)/U(k,k)
        U(j,k:n) = U(j,k:n)- U(k,k:n)*U(j,k)/U(k,k)
    end
end
D = diag(diag(U));
end

function [L,U]=directlu(A)
n = size(A); n = n(1);
L = eye(n); U = A; 
for k = 1:1:n-1
    for j = k+1:1:n
        L(j,k) = U(j,k)/U(k,k)
        U(j,k:n) = U(j,k:n)- U(k,k:n)*U(j,k)/U(k,k)
    end
end
end

function [L,U,P] = lupp(A)
n = size(A); n = n(1);
L = eye(n); U = A; P = eye(n);
for k = 1:1:n-1
    [M,i] = max(abs(U([k:n],k)))
    U([k,i+k-1],:) = U([i+k-1,k],:);
    L([k,i+k-1],1:k-1) = L([i+k-1,k],1:k-1);
    P([k,k+i-1],:) = P([k+i-1,k],:);
    for j = k+1:1:n
        L(j,k) = U(j,k)/U(k,k);
        % U(j,k) = 0;
        U(j,k:1:n) = U(j,k:1:n)- U(k,k:1:n)*U(j,k)/U(k,k);
    end
end
end

function a = reg_lin(x,fx)
n = length(x);
X = ones(n,2);
X(:,2) = x;
a = X\fx.';
end
