format long
k = 2000;
TaylorCoefficient = factorial(2*k)/(4^k*factorial(k)^2*2*k);
%% 
TaylorCoefficient = 1/(2*k);
for i = 1:k
    TaylorCoefficient = TaylorCoefficient*(k+i)/(4*i);
end
TaylorCoefficient
