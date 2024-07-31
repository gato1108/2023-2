%Load coordinates V0, C0, C1 and V1 where a given column has the (x,y)
%coordinates of a given Bezier curve. The number of columns is the number
%of Bezier curves.
load('C:\Users\Panorámica\Desktop\2023 S2\Cálculo Científico\Tarea 3\AttachedFiles\Glyph-d.mat')
t=0:0.01:1;
figure
hold on
for j=1:size(V0,2) %loop over all Bezier curves
    %Get coordinates of Bezier curve for all values in t
    v=BezierCurve(V0(:,j),C0(:,j),C1(:,j),V1(:,j),t);
    plot(v(1,:),v(2,:),'b') %plot Bezier curve
    axis equal %to visualize in scale
end

load('C:\Users\Panorámica\Desktop\2023 S2\Cálculo Científico\Tarea 3\AttachedFiles\Glyph-p.mat')
for j=1:size(V0,2) %loop over all Bezier curves
    %Get coordinates of Bezier curve for all values in t
    v=BezierCurve(V0(:,j),C0(:,j),C1(:,j),V1(:,j),t);
    v(1,:) = (v(1,:)+ones(1,101)*500);
    plot(v(1,:),v(2,:),'b') %plot Bezier curve
    axis equal %to visualize in scale
end


%Plot initial and final vertices as red asterisks
% plot(V0(1,:),V0(2,:),'r*',V1(1,:),V1(2,:),'r*')
%Plot control points as green squares.
% plot(C0(1,:),C0(2,:),'gs',C1(1,:),C1(2,:),'gs')

% load('Glyph-b.mat')

function v=BezierCurve(v0,c0,c1,v1,t)
n = length(t)-1;
v = zeros(2,n+1);
for i=1:(n+1);
    tt = t(i);
    v(:,i) = (1-tt)^3*v0+3*(1-tt)^2*tt*c0+3*(1-tt)*tt^2*c1+tt^3*v1;
end
end