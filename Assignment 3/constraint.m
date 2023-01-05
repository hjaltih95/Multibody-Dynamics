% use the symbolic toolbox to get the derivatives of the constraints
%
% by A.L. Schwab, TUDelft 2019

% define some symbolic variables

syms x1 y1 phi1 x2 y2 phi2
syms xd1 yd1 phid1 xd2 yd2 phid2 % for dx/dt I use xd etc
syms l1 l2 l0
syms omega t


% put the cm coordinates in a column vector x
x = [x1; y1; phi1; x2; y2; phi2];
% and the time derivaties
xd = [xd1; yd1; phid1; xd2; yd2; phid2];

% the constraints 
dxA = x1-l1/2*cos(phi1);
dyA = y1-l1/2*sin(phi1);
dxB = (x2-l2/2*cos(phi2))-(x1+l1/2*cos(phi1));
dyB = (y2-l2/2*sin(phi2))-(y1+l1/2*sin(phi1));

% Contact point in c assignment 2
%dC1 = x2 + l2/2*cos(phi2);

% spring Assignemnt 3
%dCl = sqrt((x1-l1*sin(phi1)/6+l1/2)^2 + (y1+l1*cos(phi1)/6)^2) - l0 

%exericse 3 Qa

% Motor places between gorund and hinge!
dCl = phi1 - omega*t;

%dCl1 = x1 + l1*cos(phi1)/2;

% vertical wall in the origin, thus the bar cannot cross this point,
% similar to point c in assignment 2 except it now depends on the position
% of the wall as well, see the 0, that is position of wall.
dCl2 = x2 + l2*cos(phi2)/2-0;
%dCl3 = y1 + l1*sin(phi1)/2 - l1;
%dCl4 = y2 + l2*sin(phi2)/2 - 2*l2;


% put all constraints in one vector C
C = [dxA; dyA; dxB; dyB; dCl2];
%Cl = [dCl4];
% Cx is the jacobian dC/dx
Cx = jacobian(C,x);
Cx = simplify(Cx)

% Cl = jacobian(Cl,x);
% Cl = simplify(Cl)

% and C2 are the convective terms C,xx*xd*xd
% which are by definition d(dC/dt)/dx*xd
% first determine dC/dt
Cd = Cx*xd % this is dC/dt=dC/dx*xd

% and next the convective terms d(dC/dt)/dx*xd
C2 = jacobian(Cd,x)*xd;
C2 = simplify(C2)



% to transform the symbolic code to matlab code either:
% copy and paste the result from teh command window in your m-file
% edit to get rid of superfluous brackets

% or use matlabFunction to make Matlab code out of it
matlabFunction(Cx,'File','Cxfunction')
matlabFunction(C2,'File','C2function')
% and then copy the code from witin the m-files into your own matlab code