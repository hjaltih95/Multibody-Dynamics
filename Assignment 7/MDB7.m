%% Set up EOM
syms alpha beta gamma
syms alphad betad gammad
syms alphadd betadd gammadd

% alpha = theta2; beta = theta4 = theta3; gamma = theta5

q=[alpha; beta; gamma];
qd=[alphad; betad; gammad];
qdd=[alphadd; betadd; gammadd];

% Given values
% Distance
O2A = 0.2;          
O4B = 0.7;          
BC = 0.6;           
O4O2 = 0.3;         
O4G4 = 0.4;         
BG5 = 0.3;          
yc = 0.9;           

% Masses
m3 = 0.5;           
m4 = 6;
m5 = 4;
m6 = 2;

% Mass moment inertia 
J4 = 10;
J5 = 6;
J2 = 100;
J3 = 0;

% Forces
F = 1000;
T = 0;

% Inital velocities
omega = 75;

% Kinematics
xA = O2A*cos(alpha);
yA = O4O2+O2A*sin(alpha);
x4 = O4G4*cos(beta);
y4 = O4G4*sin(beta);
x5 = O4B*cos(beta)+BG5*cos(gamma);
y5 = O4B*sin(beta)+BG5*sin(gamma);
x6 = O4B*cos(beta)+BC*cos(gamma);

% Transformation matrix 
Ti = [alpha; xA; yA; beta; x4; y4; beta; x5; y5; gamma; x6];   % note that theta3=theta4 so beta is used for both xA, yA and x4,y4

% Find Tij where xd = Tij*qd
Tij = jacobian(Ti,q);

% Find velocity
Tiv = Tij*qd;

% Find Tijk where xdd = Tij*qdd + Tijk*qd*qd 
Tijk=zeros(11,3);

for i=1:3
    Tijk=Tijk+jacobian(Tij(:,i),q);
end

Tacc = jacobian(Tiv,qd)*qdd + jacobian(Tiv,q)*qd;

% Convective acceleration terms
gk = Tijk*(qd.*qd);


% Mass matrix
Mij = diag([J2 m3 m3 J3 m4 m4 J4 m5 m5 J5 m6]);

% Constraints
C6 = O4G4*sin(beta) + BC*sin(gamma) - yc;
CA = xA/yA - x4/y4;

C = [C6; CA];

Cd = jacobian(C,q);
Cd = simplify(Cd);
Cdd = jacobian(Cd*qd,q)*qd;
Cdd = simplify(Cdd);

% Applied forces, no gravity and external forces
Fi = zeros(11,1);
Fi(1) = T;
Fi(11) = F;

% Find reduced mass matrix 
M = simplify(Tij.'*Mij*Tij);

Mbar = simplify([M Cd.';Cd zeros(2,2)]);

% Combined force matrix
Q = simplify(Tij.'*(Fi-Mij*gk));

Qbar = simplify([Q;-Cdd]);

Acc = Mbar\Qbar;

matlabFunction(Acc,'File','acchw7');
matlabFunction(C,'File','Chw7');
matlabFunction(Cd,'File','Cdhw7');

x6d = jacobian(x6,q)*qd;
x6dd = jacobian(x6d,qd)*qdd + jacobian(x6d,q)*qd;

matlabFunction(x6d,'File','x6dhw7');
matlabFunction(x6dd,'File','x6ddhw7');
matlabFunction(x6,'File','x6hw7');

%% Setup
time=2;
nn=13;
N=2.^nn;
h=time./N;

% Calculate initial angles and velocities
alpha = 0;
alphad = omega*2*pi/60;

beta = atan((O4O2+O2A*sin(alpha))/O2A*cos(alpha));
betad = O2A*alphad*cos(beta-alpha)/sqrt(O2A^2+O4O2^2);

gamma = pi+asin((yc-O4B*sin(beta))/BC);
%gammad = O4B*cos(beta)/(BC*sqrt(1-(yc-O4B*sin(beta))^2/BC^2));
gammad = 1.8433;

y0 = [alpha; beta; gamma; alphad; betad; gammad];
ang = [alpha; beta; gamma; alphad; betad; gammad];


% Setup for Gauss-Newton method
tol = 1e-12;
maxit = 10;
%% RK4 method with coordinate projection
    
for j=1:N
    [k1, force] = qa(y0);
    k2=qa(y0 + h/2*k1);
    k3=qa(y0 + h/2*k2);
    k4=qa(y0 + h*k3);
    qn= y0+1/6*h*(k1+2*k2+2*k3+k4);
    
    
%     % C(q_n+1)
%     tap = GNP(qn(1:3));      
%     Dc = tap(4:5);      % Set C(q_n+1)
    
%     %Find corrected positions
%     i = 0;
%     while (max(abs(Dc) > tol) || (i < maxit))
%         tap = GNP(qn(1:3));
%         dp = tap(1:3);                       % Find error
%         qn(1:3) = qn(1:3) + dp;         % Corrected position
%         tap = GNP(qn(1:3));
%         Dc = tap(4:5);                  % Next C(q_n+1)
%         i = i+1;                         %counter
%     end

    % Find corrected velocities
    tep = GNV(qn);
    dcd = tep(4:5);
    dv = tep(1:3);
    qn(4:6) = qn(4:6) + dv;
    
    y0 = qn;
    q_n(j,:)=qn;
    force_all(j,:)=force;
    
    acceleration(j,:) = [k1(4);k1(5);k1(6)];
    
    % Velocity of slider 3 and rocker 4
    vx3 = -0.2*sin(qn(1))*qn(4);
    yx3 = 0.2*cos(qn(1))*qn(4);
    v3(j,:) = sqrt(vx3^2+yx3^2);
    
    vx4 = -0.4*sin(qn(2))*qn(5);
    vy4 = 0.4*cos(qn(2))*qn(5);
    v4(j,:) = sqrt(vx4^2+vy4^2);
    v3tov4 = v3 - v4;
    
    
    
end

%% Plot the data
% b) Plot the speeds of the crank, rocker and bar
% q_plot=[ang';q_n];
% figure(1);
tt=0:h:time;
% plot(tt,q_plot(:,4)); hold on
% plot(tt,q_plot(:,5));
% plot(tt,q_plot(:,6));
% title('Speed of crank 2, rocker 4 and bar 5')
% xlabel('Time [sec]')
% ylabel('Speed rad/s')
% legend('Crank 2', 'Rocker 4', 'Bar 5')

% c)
% figure(2);
% plot(tt(2:end),v3tov4);
% title('Speed of slider 3 with respect to rocker 4')
% xlabel('Time [sec]')
% ylabel('Speed m/s')
% legend('Slider 3')



% % d)
% x6_pos = x6hw7(q_n(:,2),q_n(:,3));
% x6_vel = x6dhw7(q_n(:,5),q_n(:,2),q_n(:,6),q_n(:,3));
% x6_acc = x6ddhw7(q_n(:,5),acceleration(:,2),q_n(:,2),q_n(:,6),acceleration(:,3),q_n(:,3));
% figure(3);
% %plot(tt(2:end),x6_pos); hold on
% %plot(tt(2:end),x6_vel);
% plot(tt(2:end),x6_acc);
% title('Horizontal acceleration of slider 6')
% xlabel('Time [sec]')
% ylabel('Acceleration m/s^2')
% legend('Slider 6')


% e) and f)
figure(4)
plot(tt(2:end),force_all(:,1)); hold on
plot(tt(2:end),force_all(:,2));
title('Reaction forces of slider 6 and slider 3')
xlabel('Time [sec]')
ylabel('Force N')
legend('Slider 6', 'Slider 3')

%% Find the period of the system

ff = find(q_n(:,1)>2*pi,1,'First');
period = tt(ff);                        % 0.8780

%% Standard First-Order Form
function [acc, lambda] = qa(y)
alpha=y(1);
beta=y(2);
gamma=y(3);
alphad=y(4);
betad=y(5);
gammad=y(6);

acc = acchw7(alpha,alphad,betad,beta, gammad, gamma);

lambda = acc(4:5,1);

acc = [alphad; betad; gammad; acc(1,1); acc(2,1); acc(3,1)];


end

%% Gauus-Newton method for position coordinate projection
function delta = GNP(q)

alpha = q(1);
beta = q(2);
gamma = q(3);


D = Chw7(alpha,beta,gamma);
Dq = Cdhw7(alpha,beta,gamma);

delta = [-Dq'*inv(Dq*Dq.')*D; D];
end

%% Gauus-Newton method for velocity coordinate projection 
function delta_v = GNV(q)

alpha = q(1);
beta = q(2);
gamma = q(3);

Dq = Cdhw7(alpha,beta,gamma);

delta_v = [-Dq'*inv(Dq*Dq.')*Dq*q(4:6); Dq*q(4:6)];
end