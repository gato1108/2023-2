%% (i)

f=@(x) cos(x)
h = 1e-10:0.25e-10:1e-6;
l = length(h);
teo = (f(h+0.5)-f(0.5))./h;
error = abs(-sin(0.5)-teo);
loglog(h,error)

%% (ii)

h = 1e-3:0.5e-3:1e-1
l = length(h)
ad = (cos(h+0.5)-cos(0.5))./h;
error_ad = abs(-sin(0.5)-ad);
cen = (cos(h+0.5)-cos(0.5-h))./(2*h);
error_ce = abs(-sin(0.5)-cen);
cen5 = (-cos(2*h+0.5)+8*cos(0.5+h)-8*cos(0.5-h)+cos(0.5-2*h))./(12*h);
error_ce5 = abs(-sin(0.5)-cen5);
loglog(h,error_ad,h,error_ce,h,error_ce5)
legend('Adelantada','Centrada','5 Puntos')

%% Adelantada
p = polyfit(h_2,error_ad,6)

%% Centrada
p = polyfit(h_2,error_ce,6)

%% 5 Puntos 
p = polyfit(h_2,error_ce5,6)

%% (iii)
f=@(x) x.*abs(x)

h = 1e-3:0.5e-3:1e-1 ; 
ad = (f(h))./h;
error_ad = abs(ad);
cen = (f(h)-f(-h))./(2*h);
error_ce = abs(cen);
cen5 = (-f(2*h)+8*f(h)-8*f(-h)+f(-2*h))./(12*h);
error_ce5 = abs(cen5);
loglog(h,error_ad,h,error_ce,h,error_ce5)
legend('Adelantada','Centrada','5 Puntos')
