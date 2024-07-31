% Trapecio 
h = 1
x = -1:h:1;
f = x.^3-4*x.^2+x-1;

plot(x,f,'red')

%% Trapecio

I = 0
for j = 1:length(x)-1
    I = I+(f(j)+f(j+1))*h/2;
end
I

%% Simpson 

I = 0
for j = 1:length(x)/2
    I = I+(f(2*j-1)+f(2*j)*4+f(2*j+1))*h/3;
end
I