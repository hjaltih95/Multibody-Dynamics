%% Find holonomic and non-holonomic constraints a) and b)

syms x1 y1 phi1 x2 y2 phi2
syms x1d y1d phi1d x2d y2d phi2d
syms x1dd y1dd phi1dd x2dd y2dd phi2dd

x = [x1;y1;phi1;x2;y2;phi2];
xd = [x1d;y1d;phi1d;x2d;y2d;phi2d];
xdd = [x1dd;y1dd;phi1dd;x2dd;y2dd;phi2dd];

% Given values
a = 0.2;
b = 0.6;
c = 0.1;
d = 0.1;
m1 = 46;
m2 = 1;
I1 = 10;
I2 = 0.005;

% Position constraints in hinge B
cXB = x1 + b*cos(phi1)-x2+d*cos(phi2);
cYB = y1 + b*sin(phi1)-y2+d*sin(phi2);

Cp = [cXB;cYB];

% Velocity constraints on tires, we only look at perpendicular constraint
% for sliding, i.e. no sideslip
vA = [x1d;y1d;0] + cross([0;0;phi1d],[-a*cos(phi1);-a*sin(phi1);0]);
cVA = [-sin(phi1),cos(phi1),0]*vA;

vC = [x2d;y2d;0] + cross([0;0;phi2d],[c*cos(phi2);c*sin(phi2);0]);
cVC = [-sin(phi2),cos(phi2),0]*vC;

Cv = [cVA;cVC];

% The jacobian of the constraints
dCp = jacobian(Cp,x);
dCp = simplify(dCp);

dcV = jacobian(Cv,xd);
dcV = simplify(dcV);

Dcp = [dCp];
Dcv = [dcV];
Dc = [dCp;dcV];

% Find the convectve acceleration terms
% Holonomic constriaints
Dcpt = Dcp*xd;
ck = jacobian(Dcpt,x)*xd;
ck = simplify(ck);

% Non-holonomic constraints
Dcvt = Dcv*xd;
sk = jacobian(Dcvt,x)*xd;
sk = simplify(sk);

% conv terms
hk = [ck;sk];

matlabFunction(Dcp,'File','Dcphw8');
matlabFunction(Dcv,'File','Dcvhw8');
matlabFunction(Cv,'File','Cvhw8');
matlabFunction(Cp,'File','Cphw8');
matlabFunction(Dc,'File','Dchw8');
matlabFunction(hk,'File','hkhw8');
%% Compute DAE c)
% Define mass matrix
M=diag([m1,m1,I1,m2,m2,I2]);

% Left hand side
LHS = [M,Dc.';Dc,zeros(4,4)];

% Right hand side
F = zeros(6,1);
%F = [0;0;2*cos(pi*t);0;0;-2*cos(pi*t)];

RHS = [F;-hk];

acc = LHS\RHS;

matlabFunction(LHS,'File', 'LHShw8');
matlabFunction(RHS,'File', 'RHShw8');
matlabFunction(acc,'File', 'acchw8');


%% e) Numerical intergration, Runge-Kutta 4th order
 

% Calculate initial angles and velocities
x1 = a;
y1 = 0;
phi1 = 0;
x2 = a+b;  % a+b
y2 = d;      % d
phi2 = pi/2;

x1d = 1;    % 1
phi1d = 0;


% Calculated velocities
x2d = 0.5;   % 0.5
y2d = 0;
y1d = 0;
phi2d = 5;  % 5

% % Calculate initial vel
% eq1 = -x1d*sin(phi1)+y1d*cos(phi1)-phi1d/5;
% eq2 = -x2d*sin(phi2) + y1d*cos(phi2) + phi2d/10;
% eq3 = x1d - x2d - (3*phi1d*sin(0))/5 - (phi2d*sin(phi2))/10;
% eq4 = y1d - y2d + (3*phi1d*cos(0))/5 + (phi2d*cos(phi2))/10;
% sol = solve(eq1,eq2,eq3,eq4);
% 
% x2d = vpa(sol.x2d);
% y1d = vpa(sol.y1d);
% y2d = round(vpa(sol.y2d));   % We use round because cos and sin dont give
% exact solutions for 0, pi/2, pi, 3/2*pi, 2pi, etc.
% phi2d = vpa(sol.phi2d);

y0 = [x1;y1;phi1;x2;y2;phi2;x1d;y1d;phi1d;x2d;y2d;phi2d];
ang = [x1;y1;phi1;x2;y2;phi2;x1d;y1d;phi1d;x2d;y2d;phi2d];

% setup
time=1;
nn=13;
N=2.^nn;
h=time./N;

t= 0;
for j=1:N
    [k1, force] = qa(y0);
    k2=qa(y0 + h/2*k1);
    k3=qa(y0 + h/2*k2);
    k4=qa(y0 + h*k3);
    qn= y0+1/6*h*(k1+2*k2+2*k3+k4);
    
    t = t+h;
    
    % Correct position
    qn = projectposition(qn);
    
    % Correct velocities
     qn = projectspeed(qn);
    
    y0 = qn;
    q_n(j,:)=qn;
    force_all(j,:)=force;
    T(j) = t;
    
end

%% Plots
% e)
q_plot=[ang';q_n];
figure 
plot(q_plot(:,1),q_plot(:,2),q_plot(:,4),q_plot(:,5))
grid on
xlabel('X[m]');
ylabel('Y[M]');
legend('Body 1','Body 2');
title('Question e')



%% Standard First-Order Form
function [acc, lambda] = qa(q)

x1 = q(1);
y1 = q(2);
phi1 = q(3);
x2 = q(4);
y2 = q(5);
phi2 = q(6);

x1d = q(7);
y1d = q(8);
phi1d = q(9);
x2d = q(10);
y2d = q(11);
phi2d = q(12);

acc = acchw8(phi1,phi2,phi1d,phi2d,x1d,x2d,y1d,y2d);

lambda = acc(7:10,1);

acc = [x1d; y1d; phi1d; x2d ; y2d; phi2d; acc(1,1); acc(2,1); acc(3,1); acc(4,1); acc(5,1); acc(6,1)];


end

%% d) the two coordinate projection functions

%% Gauus-Newton method for position coordinate projection
function q_corr = projectposition(q)
% setup for while loop
tol = 1e-12;
i = 0;
maxit = 10;

x1 = q(1);
y1 = q(2);
phi1 = q(3);
x2 = q(4);
y2 = q(5);
phi2 = q(6);

x1d = q(7);
y1d = q(8);
phi1d = q(9);
x2d = q(10);
y2d = q(11);
phi2d = q(12);
dxvec = reshape (q(7:12),6,1);

q_corr=q;

s = Cphw8(phi1,phi2,x1,x2,y1,y2);

% Find corrected positions
while (tol<max(abs(s)) && i<maxit)
    Dv = Dcphw8(phi1,phi2);
    s = Cphw8(phi1,phi2,x1,x2,y1,y2);
    
    delta_p = -Dv'*inv(Dv*Dv.')*s;
    q_corr(1:6) = q_corr(1:6) + delta_p;
    
    i = i+1;
    
end

end

%% Gauus-Newton method for velocity coordinate projection 
function q_corr = projectspeed(q)

x1 = q(1);
y1 = q(2);
phi1 = q(3);
x2 = q(4);
y2 = q(5);
phi2 = q(6);

x1d = q(7);
y1d = q(8);
phi1d = q(9);
x2d = q(10);
y2d = q(11);
phi2d = q(12);
dxvec = reshape (q(7:12),6,1);

q_corr=q;

Ds = Dchw8(phi1,phi2);
s = Ds*dxvec;

% Find correceted velocities
delta_v =  -Ds.'*((Ds*Ds.')\s);
q_corr(7:12) = vpa(q_corr(7:12)+delta_v);

end

