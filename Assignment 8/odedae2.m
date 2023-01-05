function yvecd = odedae2(t,yvec)

global a b c d m1 m2 I1 I2 CD FC
% copy state vector into meaningful variables
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

A =[         46,         0,                0,          0,         0,             0,                1,               0, -sin(phi1),          0;
          0,        46,                0,          0,         0,             0,                0,               1,  cos(phi1),          0;
          0,         0,               10,          0,         0,             0, -(3*sin(phi1))/5, (3*cos(phi1))/5,       -1/5,          0;
          0,         0,                0,          1,         0,             0,               -1,               0,          0, -sin(phi2);
          0,         0,                0,          0,         1,             0,                0,              -1,          0,  cos(phi2);
          0,         0,                0,          0,         0,         1/200,    -sin(phi2)/10,    cos(phi2)/10,          0,       1/10;
          1,         0, -(3*sin(phi1))/5,         -1,         0, -sin(phi2)/10,                0,               0,          0,          0;
          0,         1,  (3*cos(phi1))/5,          0,        -1,  cos(phi2)/10,                0,               0,          0,          0;
 -sin(phi1), cos(phi1),             -1/5,          0,         0,             0,                0,               0,          0,          0;
          0,         0,                0, -sin(phi2), cos(phi2),          1/10,                0,               0,          0,          0];
 B =[ 0;
                                                0;
                                                0;
                                                2*cos(pi*t);
                                                0;
                                                -2*cos(pi*t);
 (3*cos(phi1)*dphi1^2)/5 + (cos(phi2)*dphi2^2)/10;
 (3*sin(phi1)*dphi1^2)/5 + (sin(phi2)*dphi2^2)/10;
        dphi1*dy1*sin(phi1) + dphi1*dx1*cos(phi1);
        dphi2*dy2*sin(phi2) + dphi2*dx2*cos(phi2)];
  
  DAE = A\B;
yvecd = yvec;
% copy velocties and accelerations into yd
yvecd(1:6) = yvec(7:12);
yvecd(7:12) = DAE(1:6);
% Copy constraint force in global FC
FC = DAE(4);
end
 