x1 = [.478,.781,1.069,1.438,1.914,2.92,4.05,5.14,7.15,8.9];
fx1 = [4.82,7.88,10.78,14.5,19.29,29.35,40.5,51.4,71.3,88.5];

a1 = reg_lin(fx1,x1)

yplot = [0:.01:90];
xplot = a1(1)+a1(2)*yplot;
scatter(fx1,x1,'black','filled')
hold on
plot(yplot,xplot,'red')
ylabel('Voltaje (V)')
xlabel('Corriente (mA)')
title('Resistencia')
grid on

sqrt(error(a1,fx1,x1))/10

%%

x2 = [.213,.666,1.037,1.493,2,2.54,3.24,3.82,4.26,4.91,5.11,5.84,6,6.4,7.03];
fx2 = [28.96,46,55.7,66.5,77.5,88,100.4,110,116.7,126.4,129.1,139.2,141.4,146.6,154.5];

a21 = reg_lin(fx2(1:6),x2(1:6))
a22 = reg_lin(fx2(7:15),x2(7:15))

yplot1 = [28:.01:100];
yplot2 = [88:.01:160];
xplot1 = a21(1)+a21(2)*yplot1;
xplot2 = a22(1)+a22(2)*yplot2;
scatter(fx2,x2,'black','filled')
hold on
plot(yplot1,xplot1,'red')
plot(yplot2,xplot2,'red')
ylabel('Voltaje (V)')
xlabel('Corriente (mA)')
title('Ampolleta')
grid on

sqrt((error(a21,fx2(1:6),x2(1:6))+error(a22,fx2(7:15),x2(7:15))))/15
%%

x3 = [1.506,1.868,2.25,2.64,2.84,3.21,3.38,3.56,3.71,3.92,4.14,4.33,4.45,4.78,4.91];
fx3 = [.03,1.21,3.34,5.7,6.96,9.23,10.28,11.39,12.36,13.65,15.1,16.26,17,19.14,19.95];

a31 = reg_lin(fx3(1:2),x3(1:2))
a32 = reg_lin(fx3(3:5),x3(3:5))
a33 = reg_lin(fx3(6:15),x3(6:15))

yplot1 = [0:.01:2];
yplot2 = [2:.01:9];
yplot3 = [9:.01:20];
xplot1 = a31(1)+a31(2)*yplot1;
xplot2 = a32(1)+a32(2)*yplot2;
xplot3 = a33(1)+a33(2)*yplot3;
scatter(fx3,x3,'black','filled')
hold on
plot(yplot1,xplot1,'red')
plot(yplot2,xplot2,'red')
plot(yplot3,xplot3,'red')
ylabel('Voltaje (V)')
xlabel('Corriente (mA)')
title('Diodo Led')
grid on

sqrt((error(a31,fx3(1:2),x3(1:2))+error(a32,fx3(3:5),x3(3:5))+error(a33,fx3(6:15),x3(6:15))))/15


%% LAB 3

tiempo = 0:30:600
t_sin_lana = [20,21,21,21,22,22,23,23,23,23,24,24,25,25,25,25,26,26,26,27,27]
t_con_lana = [20,22,22,23,23,23,24,24,25,25,26,26,26,27,27,27,28,28,28,29,29]

a1 = reg_lin(tiempo,t_sin_lana)
a2 = reg_lin(tiempo,t_con_lana)

xplot = [0:.01:600];
subplot(2,1,1);
yplot = a1(1)+a1(2)*xplot;
scatter(tiempo,t_sin_lana,'black','filled')
hold on
plot(xplot,yplot,'red')
ylabel('Temperatura (°C)')
xlabel('Tiempo (s)')
title('Temperatura v/s Tiempo (sin lana)')
grid on

subplot(2,1,2); 
yplot = a2(1)+a2(2)*xplot;
scatter(tiempo,t_con_lana,'black','filled')
hold on
plot(xplot,yplot,'red')
ylabel('Temperatura (°C)')
xlabel('Tiempo (s)')
title('Temperatura v/s Tiempo (con lana)')
grid on

function a = reg_lin(x,fx)
n = length(x);
X = ones(n,2);
X(:,2) = x;
a = X\fx.';
end

function e = error(a,x,fx)
n = length(fx);
e = 0;
for i=1:n
    e = e+(a(1)+a(2)*x(i)-fx(i))^2;
end
end