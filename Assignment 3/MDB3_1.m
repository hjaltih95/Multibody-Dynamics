%% Find derivatives of constraints a)
% Set up variables
syms x1 y1 phi1 x2 y2 phi2
syms xd1 yd1 phid1 xd2 yd2 phid2 % for dx/dt I use xd etc
syms l1 l2 l0 ls omega t

% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2]
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2]

% the constraints 
dxA = x1-l1/2*cos(phi1);
dyA = y1-l1/2*sin(phi1);
dxB = (x2-l2/2*cos(phi2))-(x1+l1/2*cos(phi1));
dyB = (y2-l2/2*sin(phi2))-(y1+l1/2*sin(phi1));

dC1 = phi1-omega*t
dC2=x2+l2/2*cos(phi2)
% put all constraints in one vector C
C = [dxA; dyA; dxB; dyB; dC1];
Cs = sqrt((x1+1/6*l1*cos(phi1)+l1/2)^2+(y1+1/6*l1*sin(phi1)-0)^2)-l0;

% Cx is the jacobian dC/dx
Cx = jacobian(C,x);
Cx = simplify(Cx)

% Csx is for the spring constraint
Csx = jacobian(Cs,x);
Csx = simplify(Csx)

% and C2 are the convective terms C,xx*xd*xd
% which are by definition d(dC/dt)/dx*xd
% first determine dC/dt
Cd = Cx*xd; % this is dC/dt=dC/dx*xd

% and next the convective terms d(dC/dt)/dx*xd
C2 = jacobian(Cd,x)*xd;
C2 = simplify(C2)

matlabFunction(Cx,'File','Cxfunction')
matlabFunction(C2,'File','C2function')