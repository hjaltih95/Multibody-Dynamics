function [conv,sc,Ds,h] = constraint(yvec)

global A_ M1_ M2_ I1_ I2_ F_ FC_ CD_ B_ C_ D_

% copy state vector into meaningful variables

x1 = yvec(1);
y1 = yvec(2);
p1 = yvec(3);
x2 = yvec(4);
y2 = yvec(5);
p2 = yvec(6);
dx1 = yvec(7);
dy1 = yvec(8);
dp1 = yvec(9);
dx2 = yvec(10);
dy2 = yvec(11);
dp2 = yvec(12);
dxvec = reshape(yvec(7:12),6,1);

%Constraint Matrix

Ds=[-sin(p1), cos(p1), -A_ , 0, 0, 0;
0, 0, 0, -sin(p2), cos(p2), C_ ;
1, 0, -B_ *sin(p1), -1, 0, -D_*sin(p2);
0, 1, B_ *cos(p1), 0, -1, D_*cos(p2)];

% Holonomic constraints

pos = [x1+B_*cos(p1)+D_*cos(p2)-x2;
y1+B_*sin(p1)+D_*sin(p2)-y2];

% Nonholonomic constraints

s = [Ds(1,:);Ds(2,:)]*dxvec;

% Convective terms

conv=[-B_*cos(p1)*dp1^2 - D_*cos(p2)*dp2^2;
 -B_*sin(p1)*dp1^2 - D_*sin(p2)*dp2^2];

%Nonholonomic terms

h = [  -dp1*dx1*cos(p1) - dp1*dy1*sin(p1);
-dp2*dx2*cos(p2) - dp2*dy2*sin(p2)];
sc=[s;pos];
end