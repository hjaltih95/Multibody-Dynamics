
%% Set up EoM a)
% Set up variables
syms alpha beta gamma
syms alphadot betadot % for dx/dt I use xd etc
syms alpha2dot beta2dot
syms l m M I g Qa Qb

q=[alpha; beta];
qd=[alphadot; betadot];

% Define CoM coordinates
x1=sin(alpha+gamma)*l/2;        % x coordinates of stance leg
y1=cos(alpha+gamma)*l/2;        % y coordinates of stance leg
phi1=alpha+gamma;               % angle of stance leg
xh=sin(alpha+gamma)*l;          % x coordinates of hips
yh=cos(alpha+gamma)*l;          % x coordinates of hips
x2=xh + sin(beta-gamma)*l/2;    % x coordinates of swing leg
y2=yh - cos(beta-gamma)*l/2;    % x coordinates of swing leg
phi2=beta-gamma;                % angle of swing leg


% Coordinate matrix
X=[x1;y1;phi1;xh;yh;x2;y2;phi2];
Xij=simplify(jacobian(X,q));
% velocity
Xd=Xij*qd;

Tijk=zeros(8,2);
for i=1:2
    Tijk=Tijk+jacobian(Xij(:,i),q);
end
% find acceleration
xdd = simplify(Xij*[alpha2dot;beta2dot] + Tijk*(qd.*qd));

% Define mass matrix
Mij=diag([m, m, I, M, M, m, m, I]);

% Define potential energy
V=-(m*g*y1 + m*g*y2+ M*g*yh);

% Define kinetic energy
T = simplify(1/2*Xd.'*Mij*Xd);

% External forces
Q = [Qa;Qb];

% Find lagrange 
Vdq = simplify(jacobian(V,q));
Tdq = simplify(jacobian(T,q));
Tdqd = simplify(jacobian(T,qd));
Tdqdd = simplify(jacobian(Tdqd,q));


% Right hand side of the matrix
F_l = simplify(Q-Vdq.' + Tdq.' - Tdqdd*qd);

% Left hand side of the matrix
M_l = jacobian(Tdqd,qd);


%% Set up EoM b)
% Set up variables
syms alpha beta gamma
syms alphadot betadot % for dx/dt I use xd etc
syms alpha2dot beta2dot
syms l m M I g Qa Qb

q=[alpha; beta];
qd=[alphadot; betadot];

% Define CoM coordinates
x1=sin(alpha+gamma)*l/2;        % x coordinates of stance leg
y1=cos(alpha+gamma)*l/2;        % y coordinates of stance leg
phi1=alpha+gamma;               % angle of stance leg
xh=sin(alpha+gamma)*l;          % x coordinates of hips
yh=cos(alpha+gamma)*l;          % x coordinates of hips
x2=xh + sin(beta-gamma)*l/2;    % x coordinates of swing leg
y2=yh - cos(beta-gamma)*l/2;    % x coordinates of swing leg
phi2=beta-gamma;                % angle of swing leg

% Define transformation matrix
Ti = [x1; y1; phi1; xh; yh; x2; y2; phi2];

% Find Tij, where xd = Tij*qd
Tij = jacobian(Ti,q);

% Find Tijk where xdd = Tij*qdd + Tijk*qd*qd 
Tijk=zeros(8,2);

for i=1:2
    Tijk=Tijk+jacobian(Tij(:,i),q);
end
% Mass matrix
Mij=diag([m, m, I, M, M, m, m, I]);

% applied forces
Fi=Mij*[0,g,0,0,g,0,g,0].';

% External forces
Q = [Qa;Qb];

% Find reduced mass matrix
Mbar = simplify(Tij.'*Mij*Tij);

% Convective acceleration terms
gk=Tijk*(qd.*qd);

% Combined force matrix
Qbar = simplify(Q+Tij.'*(Fi-Mij*gk));

