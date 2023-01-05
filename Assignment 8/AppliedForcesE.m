function f = AppliedForcesE(t,yvec)
global a b c d m1 m2 I1 I2 CD FC F M
% copy state vector into meaningful variables
x1 = yvec(1);
y1 = yvec(2);
phi1 = yvec(3);
x2 = yvec(4);
y2 = yvec(5);
phi2 = yvec(6);
x1d = yvec(7);
y1d = yvec(8);
phi1d = yvec(9);
x2d = yvec(10);
y2d = yvec(11);
phi2d = yvec(12);

M0 = 2;
omega = pi;
M1 = M0 * cos(omega *t);
f = zeros(6,1);
F = [0, 0, 0, 0, 0, 0].';  % e)
% add torques questions: f)-j)
f(1) = F(1);
f(2) = F(2);
f(3) = M1;
f(4) = F(4);
f(5) = F(5);
f(6) = -M1;
end