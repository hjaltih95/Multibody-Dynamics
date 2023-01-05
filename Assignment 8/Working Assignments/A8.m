clear all
clc

syms a b c d real
syms x1 y1 x2 y2
syms x1d y1d x2d y2d
syms x1dd y1dd x2dd y2dd
syms phi1 phi1d phi1dd phi2 phi2d phi2dd
syms m1 I1 m2 I2

global a b c d m1 I1 m2 I2

x = [x1 ;y1; phi1; x2; y2; phi2];
xd = [x1d ;y1d; phi1d; x2d; y2d; phi2d];
xdd = [x1dd; y1dd; phi1dd;x2dd; y2dd; phi2dd];

% Param values




%% Holonomic COnstraints

Cx = x1 + b*cos(phi1) - x2 +d*cos(phi2);
Cy = y1 +b*sin(phi1) - y2 +d*sin(phi2);

Cm = [Cx;Cy];
Cmj = jacobian(Cm,x);
Cmi = Cmj.';
Cmjl_qdqd = jacobian(Cmj*xd,x)*xd;

%% Non Holonimoc Constraints

S1Aj = [-sin(phi1); cos(phi1);-a;0;0;0];
S1Cj = [0;0;0;-sin(phi2);cos(phi2);c];
Smi = [S1Aj S1Cj];
Smj = Smi.';

Skjl_dqdq = jacobian(Smj*xd,x)*xd

%% Forces
Fi = [zeros(6,1)];

%% Mass 
Mij = diag([m1 m1 I1 m2 m2 I2]);

%% EOM

X = [Mij Cmi Smi; Cmj zeros(2,4) ; Smj zeros(2,4)];
Y = [Fi; -Cmjl_qdqd;-Skjl_dqdq];
%sol = simplify(X\Y)
