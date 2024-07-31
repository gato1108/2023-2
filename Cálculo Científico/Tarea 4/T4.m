n = 100
x = -1:2/n:1 
% x = cos((2*(0:n)+1)*pi/(2*(n+1)))
% x = cos((0:n)*pi/n)

cond(Van(n+1,x))
cond(Che(n+1,x))
%%
n = 500
x = cos((2*(0:n)+1)*pi/(2*(n+1)));
fx_1 = abs(x);
C = Che(n+1,x);
c_1 = C\fx_1.';
fx_2 = abs(x).^7;
c_2 = C\fx_2.';
fx_3 = (1+25*x.^2).^(-1);
c_3 = C\fx_3.';

x_log = 0:1:n
semilogx(x_log,abs(c_1),x_log,abs(c_2),x_log,abs(c_3))
legend('|x|','|x|^7','1/(1+25x^2)')
grid on
% num = [0:1:n] 
% dot(chebyshevT(num, 1),c)

%%
 x = -1:0.001:1;
 err = sin(x)-(x-(1/6)*x.^3);
 max(abs(err))
 1/120

%%
 x = -1:0.001:1;
 err = sin(x)-(x*(383/384)-(5/32)*x.^3);
 max(abs(err))
1/(120*16)+1/factorial(7)
 
%%

n = 20
x = cos((2*(0:n)+1)*pi/(2*(n+1)));
fx_1 = abs(x);
C = Che(n+1,x);
c_1 = C\fx_1.';
fx_2 = abs(x).^7;
c_2 = C\fx_2.';
fx_3 = (1+25*x.^2).^(-1);
c_3 = C\fx_3.';


num = [0:1:n] 
x = -1:0.01:1;
fx = zeros(201,1);
for i=1:201
    fx(i)=dot(chebyshevT(num, x(i)),c_2);
end

plot(x,fx)
%%

function V = Van(n,x)
V = zeros(n,n);
for i=1:n
    for j=1:n
        V(i,j)=x(i)^(j-1);
end
end
end

function V = Che(n,x)
V = zeros(n,n);
for i=1:n
    for j=1:n
        V(i,j)=cos((j-1)*acos(x(i)));
end
end
end