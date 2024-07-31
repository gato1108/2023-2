%% (B)
n = 10;
xnod = -1:(2/n):1
cheb = ones(1,n+1)
for i = 1:(n+1)
    cheb(i) = -cos(pi*(i-1)/n)
end
w_xnod = BarycentricWeights(xnod);
w_cheb = BarycentricWeights(cheb);
plot(xnod,w_xnod,'blue',cheb,w_cheb,'red')
legend('Equidistantes','Chebyshev')
% print -deps epsFig

%% C

x_ev = -1:1/100:1;
y_ev_eq = ones(1,11);
y_ev_ch = ones(1,11);
y_ev = ones(1,201)
for i = 1:11
    y_ev_eq(i) = fun(xnod(i)) ;
    y_ev_ch(i) = fun(cheb(i)) ;
end

for i = 1:201
    y_ev(i) = fun(x_ev(i)) ;
end

equidistantes = BarycentricEvaluation(xnod,w_xnod,y_ev_eq,x_ev);
chebyshev = BarycentricEvaluation(cheb,w_cheb,y_ev_ch,x_ev);

error_eq = abs(y_ev-equidistantes);
error_ch = abs(y_ev-chebyshev);

plot(x_ev,error_eq,'blue',x_ev,error_ch,'red')
legend('Error Equidistantes','Error Chebyshev')
ylabel('Error') 

%% D
n = 10;
xnod = -1:(2/n):1
cheb = ones(1,n+1)
for i = 1:(n+1)
    cheb(i) = -cos(pi*(i-1)/n)
end
w_xnod = BarycentricWeights(xnod);
w_cheb = BarycentricWeights(cheb);

x_ev = -1:1/100:1;
y_ev_eq = ones(1,11);
y_ev_ch = ones(1,11);
y_ev = ones(1,201)
for i = 1:11
    y_ev_eq(i) = abs(xnod(i)) ;
    y_ev_ch(i) = abs(cheb(i)) ;
end

for i = 1:201
    y_ev(i) = abs(x_ev(i)) ;
end

equidistantes = BarycentricEvaluation(xnod,w_xnod,y_ev_eq,x_ev);
chebyshev = BarycentricEvaluation(cheb,w_cheb,y_ev_ch,x_ev);

error_eq = abs(y_ev-equidistantes);
error_ch = abs(y_ev-chebyshev);

plot(x_ev,error_eq,'blue',x_ev,error_ch,'red')
legend('Error Equidistantes','Error Chebyshev')
ylabel('Error') 

%% FUNCIONES
function w = BarycentricWeights(xnod)
n = length(xnod)
w = ones(1,n)
for i = 1:n
    prod = 1;
    for j = 1:n
        if j~=i
            prod = prod/(xnod(i)-xnod(j));
        end
    end
    w(i) = prod;
end
end

function P = BarycentricEvaluation(xnod,w,f,xeval)
n = length(xnod);
m = length(xeval);
P = ones(1,m);
for i=1:m
    x = xeval(i);
    if ismember(x, xnod)
        for j = 1:n
            if x==xnod(j)
                P(i)=f(j);
            end
        end
    else 
        num = 0;
        den = 0;
        for j = 1:n
            num = num+w(j)*f(j)/(x-xnod(j));
            den = den+w(j)/(x-xnod(j));
        end
        P(i)=num/den;
    end
end
end

function res = fun(x)
res = 1/(1+25*x^2);
end