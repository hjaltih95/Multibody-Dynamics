function dyvec = odedae(t,yvec)

global a b c d m1 m2 I1 I2 M FC

x1 = yvec(1); 
y1 = yvec(2); 
phi1 = yvec(3);
x2 = yvec(4); 
y2 = yvec(5); 
phi2 = yvec(6);

dx1 = yvec(7); 
dy1 = yvec(8); 
dphi1 = yvec(9); 
dx2 = yvec(10); 
dy2 = yvec(11); 
dphi2 = yvec(12);

M = diag([m1 m1 I1 m2 m2 I2]);
[s,Ds,h] = constraint(yvec);
DAE = [M Ds.';
       Ds zeros(4)];

% determine the applied forces
f = AppliedForces(t,yvec);
rhs = [f; -h];
lhs = DAE\rhs;
dyvec = yvec;

% copy velocties and accelerations into dy
dyvec(1:6) = yvec(7:12);
dyvec(7:12) = lhs(1:6);
% copy constraint force in global FC
FC = lhs(7:10);
end

