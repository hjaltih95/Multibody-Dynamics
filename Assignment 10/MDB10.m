%% CoM positions

syms alpha beta gamma dalpha dbeta dgamma ddalpha ddbeta ddgamma
% syms e d m1 m2 g

d = 0.3;
e = 0.375;

m1 = 1.9;
m2 = 1.5;

g = 9.81;

% Rotation matrices
R_alpha = [cos(alpha) 0 sin(alpha);...
            0 1 0;...
            -sin(alpha) 0 cos(alpha)];
        
R_beta = [1 0 0;...
          0 cos(beta) -sin(beta);...
          0 sin(beta) cos(beta)];
      
R_gamma = [cos(gamma) 0 sin(gamma);...
           0 1 0;...
           -sin(gamma) 0 cos(gamma)];
       
% Upper arm
x_u = R_alpha*R_beta*[0;0;-d/2];

% Forarm
x_f = R_alpha*R_beta*[0;0;-d] + R_alpha*R_beta*R_gamma*[4/10*e;0;0];

% Transformation matrix
Ti = [x_u;alpha;beta;x_f;gamma];

T1 = [x_u;alpha;beta];

T2 = [x_f;gamma];


%% Angular velocities

Bomega_u = [dbeta;0;0] + R_beta.'*[0;dalpha;0];

% This gives matrix
B_u = [0 1 0;...
       cos(beta) 0 0;...
       -sin(beta) 0 0];


Bomega_f = [0;dgamma;0] + R_gamma.'*[dbeta;0;0] + R_gamma.'*R_beta.'*[0;dalpha;0];

% This gives matrix 
B_f = [sin(beta)*sin(gamma) cos(gamma) 0;...
        cos(beta) 0 1;...
        -cos(gamma)*sin(beta) sin(gamma) 0];
    

%% EoM f)
% Mass moment of inertia
% Upper arm
Ixx_u = 0.015;
Iyy_u = 0.014;
Izz_u = 0.002;

% Forarm
Ixx_f = 0.001;
Iyy_f = 0.019;
Izz_f = 0.019;

Ib = [Iyy_u 0 0;...
      0 Ixx_u 0;...
      0 0 Iyy_f];
  
Ib1 = [Ixx_u 0 0;...
       0 Iyy_u 0;...
       0 0 Izz_u];
   
Ib2 = [Ixx_f 0 0;...
       0 Iyy_f 0;...
       0 0 Izz_f];

% mass matrix
M = diag([m1 m1 m1 0 0 m2 m2 m2 0]);

M1 = diag([m1 m1 m1 0 0]);

M2 = diag([m2 m2 m2 0]);


% Make T
q = [alpha;beta;gamma];
qd = [dalpha;dbeta;dgamma];
qdd = [ddalpha;ddbeta;ddgamma];

% Find Tij where xd = Tij*qd
Tij = jacobian(Ti,q);
Tij = simplify(Tij);
matlabFunction(Tij,'File','Tijhw10');

Tij1 = jacobian(T1,q);
Tij1 = simplify(Tij1);

Tij2 = jacobian(T2,q);
Tij2 = simplify(Tij2);

% Find velocity
Tiv = Tij*qd;
matlabFunction(Tiv,'File','Tivhw10');

% Find convective terms
h = jacobian((Tij*qd),q)*qd;
matlabFunction(h,'File','hhw10');

h1 = jacobian((Tij1*qd),q)*qd;

h2 = jacobian((Tij2*qd),q)*qd;


% Find acceleration, where xdd = Tij*qdd + h
Tacc = Tij*qdd + h;
matlabFunction(Tacc,'File','Tacchw10');

% Find reduced mass matrix, Tij'*M*Tij
Mbar = Tij1.'*M1*Tij1 + Tij2.'*M2*Tij2 + B_u.'*Ib1*B_u + B_f.'*Ib2*B_f;

% Mbar = Tij.'*M*Tij + B_f.'*Ib*B_f;
Mbar = simplify(Mbar);
matlabFunction(Mbar,'File','Mbarhw10');

% Find applied forces
Fi = M*[0;0;-g;0;0;0;0;-g;0]; 
Fi1 = M1*[0;0;-g;0;0];
Fi2 = M2*[0;0;-g;0];

% Torques applied to keep stabalized
% angles = [0;0;0;0;0;0];
angles = [deg2rad(30);deg2rad(20);deg2rad(60);0;0;0];   %Change angles as needed for torque
Q = torq(angles);
% Q=0;

% applied torques
Mt = 0;

% Convective terms for Euler
g1 = jacobian(B_u*qd,q)*qd;
g2 = jacobian(B_f*qd,q)*qd;


% Find reduced force matrix
Qbar = Tij1.'*(Fi1-M1*h1) + Tij2.'*(Fi2-M2*h2) + B_u.'*(Mt-Ib1*g1 - cross(Bomega_u,(Ib1*Bomega_u))) + B_f.'*(Mt-Ib2*g2 - cross(Bomega_f,(Ib2*Bomega_f))) + Q;
Qbar = simplify(Qbar);
matlabFunction(Qbar,'File','Qbarhw10');
%%
% Acceleration 
acc = Mbar\Qbar;

matlabFunction(acc,'File','acchw10');

%% Find torque
% Torgue needed for initial position, all angles and angular velocities
% zero
ang = [0;0;0;0;0;0];

torques = torq(ang);

%% A few configurations
% Configuration 1, initial position
a1 = [0;0;0;0;0;0];
alpha = a1(1);
beta = a1(2);
gamma = a1(3);
dalpha = a1(4);
dbeta = a1(5);
dgamma = a1(6);
acc1 = acchw10(alpha,beta,dalpha,dbeta,dgamma,gamma);
ddalpha = acc1(1);
ddbeta = acc1(2);
ddgamma = acc1(3);

accCoM1 = Tacchw10(alpha,beta,dalpha,dbeta,ddalpha,ddbeta,ddgamma,dgamma,gamma);


% Configuration 2, arm down
a2 = [0;0;pi/2;0;0;0];
alpha = a2(1);
beta = a2(2);
gamma = a2(3);
dalpha = a2(4);
dbeta = a2(5);
dgamma = a2(6);
acc2 = acchw10(alpha,beta,dalpha,dbeta,dgamma,gamma);
ddalpha = acc2(1);
ddbeta = acc2(2);
ddgamma = acc2(3);

accCoM2 = Tacchw10(alpha,beta,dalpha,dbeta,ddalpha,ddbeta,ddgamma,dgamma,gamma);

% Configuration 3, arm up
a3 = [pi;0;pi/2;0;0;0];
alpha = a3(1);
beta = a3(2);
gamma = a3(3);
dalpha = a3(4);
dbeta = a3(5);
dgamma = a3(6);
acc3 = acchw10(alpha,beta,dalpha,dbeta,dgamma,gamma);
ddalpha = acc3(1);
ddbeta = acc3(2);
ddgamma = acc3(3);

accCoM3 = Tacchw10(alpha,beta,dalpha,dbeta,ddalpha,ddbeta,ddgamma,dgamma,gamma);

%% g) Reduced mass matrix for all angles zero
alpha = 0;
beta = 0;
gamma = 0;

Mbar = Mbarhw10(beta,gamma);

%% h)

alpha = deg2rad(-70);
beta = deg2rad(70);
gamma = deg2rad(-30);

dalpha = 0;
dbeta = 0;
dgamma = 0;

ang = [alpha;beta;gamma;dalpha;dbeta;dgamma];

torques_ball = torq(ang);


%% i) numerical integration, Runge-Kutte 4th order

% setup
time = 5;
nn=6;
N=2^nn;
h = time./N;

t = 0;

% Initialize angles and angular velocities
alpha = deg2rad(30);
beta = deg2rad(20);
gamma = deg2rad(60);

dalpha = 0;
dbeta = 0;
dgamma = 0;

y0 = [alpha;beta;gamma;dalpha;dbeta;dgamma];
ang = [alpha;beta;gamma;dalpha;dbeta;dgamma];


for i=1:N
    k1 = state(y0);
    k2 = state(y0+h/2*k1);
    k3 = state(y0+h/2*k2);
    k4 = state(y0+h*k3);
    qn = y0 + 1/6*h*(k1+2*k2+2*k3+k4);
    
    y0 = qn;
    q_n(i,:) = qn;
    t = t+h;
    T(i) = t;
end




%% Plots
q_plot = [ang';q_n];
T_plot = [0;T'];

figure(1)
plot(T_plot,q_plot(:,1),T_plot,q_plot(:,2),T_plot,q_plot(:,3))
xlabel('Time [s]')
ylabel('Euler angles [rad]')
title('Euler angles as a function of time')
legend('\alpha','\beta','\gamma')


figure(2)
plot(T_plot,q_plot(:,4),T_plot,q_plot(:,5),T_plot,q_plot(:,6))
xlabel('Time [s]')
ylabel('ANgular velocities [rad/s]')
title('angular velocities as a function of time')
legend('\omega_\alpha','\omega_\beta','\omega_\gamma')






%% Functions

% Finding torque needed for static equlibrium
function torques = torq(ang)
% Find torques needed
alpha = ang(1);
beta = ang(2);
gamma = ang(3);
dalpha = ang(4);
dbeta = ang(5);
dgamma = ang(6);

% setup, parameters
m1 = 1.9;
m2 = 1.5;

g = 9.81;

% Mass matrix
M = diag([m1 m1 m1 0 0 m2 m2 m2 0]);

% Applied forces
Fi = -M*[0;0;-g;0;0;0;0;-g;0]; 

% convective terms
h = hhw10(alpha,beta,dalpha,dbeta,dgamma,gamma);
Tij = Tijhw10(alpha,beta,gamma);

% Torques equal the right hand side
torques = Tij.'*(Fi-M*h);
end

function states = state(y)
alpha = y(1);
beta = y(2);
gamma = y(3);
dalpha = y(4);
dbeta = y(5);
dgamma = y(6);

states = acchw10(alpha,beta,dalpha,dbeta,dgamma,gamma);

states = [dalpha;dbeta;dgamma;states(1);states(2);states(3)];

end
