function [s,Ds, h] = constraint(yvec)
% Holonomic and non-holonomic constraints s, Jacobian Ds, convective term h
global a b c d m1 m2 I1 I2 FC 

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

dxvec = reshape(yvec(7:12),[6,1]);

Ds = [1, 0, -b*sin(phi1), -1, 0, -d*sin(phi2);
      0, 1, b*cos(phi1), 0, -1, d*cos(phi2);
      -sin(phi1), cos(phi1), -a, 0, 0, 0;
      0, 0, 0, -sin(phi2), cos(phi2), c
      ];
  
 % constraint equations
 s = Ds*dxvec;
 
 % convective terms for holonomic and non-holonomic
h = [-b*cos(phi1)*dphi1^2 - d*cos(phi2)*dphi2^2;
-b*sin(phi1)*dphi1^2 - d*sin(phi1)*dphi2^2;
-dphi1*(dx1*cos(phi1) + dy1*sin(phi1));
-dphi2*(dx2*cos(phi2) + dy2*sin(phi2))];

end

