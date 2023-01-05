function dyvec = odedae2(t,yvec)
% DAE for EzyRoller
% wb1413, Multibody Dynamics B, Spring 2016.
% Arend L. Schwab 23-Apr-2016
% Copyright (c) 2016 by TU-Delft, the Netherlands.
% Modified by Iris Korevaar #4085205
global a b c d m1 m2 I1 I2 CD FC
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
% DAE
A =[ 1, 0, 0, 0, 0, 0, 1,0, -sin(phi1), 0;
0, 1, 0, 0, 0, 0, 0,1, cos(phi1), 0;
0, 0, 1/10, 0, 0, 0, -sin(phi1)/2, cos(phi1)/2,-1/2, 0;
0, 0, 0, 0, 0, 0, -1,0, 0, -sin(phi2);
0, 0, 0, 0, 0, 0, 0,-1, 0, cos(phi2);
0, 0, 0, 0, 0, 0, -sin(phi2)/8, cos(phi2)/8,0, 1/8;
1, 0, -sin(phi1)/2, -1, 0, -sin(phi2)/8, 0,0, 0, 0;
0, 1, cos(phi1)/2, 0, -1, cos(phi2)/8, 0,0, 0, 0;
-sin(phi1), cos(phi1), -1/2, 0, 0, 0, 0,0, 0, 0;
0, 0, 0, -sin(phi2), cos(phi2), 1/8, 0,0, 0, 0];
B = [ 0;
0;
2*cos(pi*t);
0;
0;
-2*cos(pi*t);
(cos(phi1)*phi1d^2)/2 + (cos(phi2)*phi2d^2)/8;
(sin(phi1)*phi1d^2)/2 + (sin(phi2)*phi2d^2)/8;
phi1d*(x1d*cos(phi1) + y1d*sin(phi1));
phi2d*(x2d*cos(phi2) + y2d*sin(phi2))];
DAE = A\B;
dyvec = yvec;
% copy velocties and accelerations into yd
dyvec(1:6) = yvec(7:12);
dyvec(7:12) = DAE(1:6);
% Copy constraint force in global FC
FC = DAE(4);
end
