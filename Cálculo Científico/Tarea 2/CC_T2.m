format long

digits(100000)

pi_str = char(vpa(pi));
ana = pi_str(2+23:2+32)

% e
% _str = char(vpa(exp(1)));
% dav_1 = e_str(2+51226:2+51230)

r_str = char(vpa(sqrt(2)));
dav_2 = r_str(2+23859:2+23854)

dav = dav_1+dav_2

%% 1

%% D
x = vpa(exp(1),100000);
disp(x);
%% 4

xnod = [1,2,7,9]
w = BarycentricWeights(xnod)

%% 
function P = BarycentricEvaluation(xnod,w,f,xeval)
P
end

%%

BarycentricEvaluation(xnod)

function w = BarycentricWeights(xnod)
n = length(xnod)
w = ones(1,n)
for i = 1:n
    prod = 1
    for j = 1:n
        if j~=i
            prod = prod/(xnod(i)-xnod(j))
        end
    end
    w(i) = prod
end
end

function P = BarycentricEvaluation(xnod)
P = BarycentricWeights(xnod)
P = P+1
end

%%
5+5


