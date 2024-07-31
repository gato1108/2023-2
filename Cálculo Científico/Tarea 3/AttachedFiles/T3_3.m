%% 3A
t_=0:0.01:1
v0_=[0;0];
v1_=[1;0];
c0_=[2;1];
c1_=[-1;1]

cor = BezierCurve(v0_,c0_,c1_,v1_,t_)

plot(cor(1,:),cor(2,:))

%% 3B

R = 1
H = 4/3*(sqrt(2)-1)
t_=0:0.01:1;
v0_=[R;0];
v1_=[0;R];
c0_=[R;H];
c1_=[H;R]

cor = BezierCurve(v0_,c0_,c1_,v1_,t_);

t=linspace(0,1);
circ_x = cos(t*pi/2);
circ_y = sin(t*pi/2);

plot(cor(1,:),cor(2,:),'red',circ_x,circ_y,'blue')

mas_lejano = (cor(:,1)+cor(:,2))/2
error = abs(1-norm(mas_lejano))

%%

Glyph-a.mat
%% FUNCIONES
function v=BezierCurve(v0,c0,c1,v1,t)
n = length(t)-1
v = zeros(2,n+1)
for i=1:(n+1)
    tt = t(i);
    v(:,i) = (1-tt)^3*v0+3*(1-tt)^2*tt*c0+3*(1-tt)*tt^2*c1+tt^3*v1;
end
end