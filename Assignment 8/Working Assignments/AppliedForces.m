function f = AppliedForces(t,yvec)

global A_ M1_ M2_ I1_ I2_ F_ FC_ CD_ B_ C_ D_

M0 = 2; % [Nm]
w = pi; % [rad/s]
M1 = M0*cos(w*t); % [Nm]

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

f = zeros(6,1);

% add drag forces to constant forces F

f(1) = F_(1);
f(2) = F_(2);
f(3) = M1;
f(4) = F_(4);
f(5) = F_(5);
f(6) = -M1;
end