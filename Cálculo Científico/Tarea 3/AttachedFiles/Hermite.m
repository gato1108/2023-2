function H=Hermite(xk,f,df,x)  
n=length(xk)-1; 
m=length(x);   
Lk=ones(n+1,m);
a=zeros(n+1,m);
b=zeros(n+1,m); 
H=zeros(1,m);  
for k=1:n+1 
   dL=0;       
   for j=1:n+1 
       if j~=k
           dL=dL+1/(xk(k)-xk(j));
           Lk(k,:)=Lk(k,:).*(x-xk(j))/(xk(k)-xk(j)); 
       end
   end
   Lsq=Lk(k,:).^2;
   b(k,:)=(x-xk(k)).*Lsq;          
   a(k,:)=(1-2*(x-xk(k))*dL).*Lsq;  
   H=H+a(k,:)*f(k)+b(k,:)*df(k);    
end