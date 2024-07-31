format long

%%
N = 100000;
a = zeros(1,N+1);
b = zeros(1,N+1);
c_1 = zeros(1,N+1);
c_2 = zeros(1,N+1);

for k = 0:N
    if k==0
        a(1) = 4;
        b(1) = 6/sqrt(3);
        c_1(1) = 1/5;
        c_2(1) = 1/239;
    else
        a(k+1) = (a(k)+(4*(-1)^k/(2*k+1)));
        b(k+1) = (b(k)+6*((-1)^k)*(1/sqrt(3))^(2*k+1)/(2*k+1));
        c_1(k+1) = (c_1(k)+((-1)^k)*(1/5)^(2*k+1)/(2*k+1));
        c_2(k+1) = (c_2(k)+((-1)^k)*(1/239)^(2*k+1)/(2*k+1));
    end
end
c = 16*c_1-4*c_2;
a_ = a(N+1)
b_ = b(N+1)
c_ = c(N+1)
%% (B)
a_rel = zeros(1,N+1);
b_rel = zeros(1,N+1);
c_rel = zeros(1,N+1);

for k = 1:N+1
    a_rel(k) = abs(pi-a(k))/(pi);
    b_rel(k) = abs(pi-b(k))/(pi);
    c_rel(k) = abs(pi-c(k))/(pi);
end
%% (C)
N_rango = 10:1:N+1;
loglog(N_rango,a_rel(N_rango),'black',N_rango,b_rel(N_rango),'red',N_rango,c_rel(N_rango),'blue')
xlabel('N')
ylabel('error relativo')
legend('(a)','(b)','(c)')
grid on
%% (E)
a_rel_log = log(a_rel(N_rango));
N_rango_log = log(N_rango);
p = polyfit(N_rango_log,a_rel_log,1);
pendiente = p(1)