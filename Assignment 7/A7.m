%% initial conditions and parameters

 O2A = 0.2; %[m]
 O4B = 0.7; %[m]
 BC = 0.6; %[m]
 O4O2 = 0.3; %[m]
 O4G4 = 0.4; %[m]
 BG5 = 0.3; %[m]
 yC = 0.9; %[m]

lengths = [O2A, O4B, BC, O4O2, O4G4, BG5, yC];

m3 = 0.5; %[kg]
 m4 = 6; %[kg]
 m5 = 4; %[kg]
 m6 = 2; %[kg]
 J2 = 100; %[kg m�2]
 J4 = 10; %[kg m�2]
 J5 = 6; %[kg m�2]

 % Force and Torque
 F = 1000; %[N]
 T = 0; %[Nm]
 omega_2 = 75; %[rpm CCW]

 % Initial conditions
 sub_theta2 = 0; % [rad]
sub_dtheta2_rpm = 75; % [rpm CCW]
 sub_dtheta2 = sub_dtheta2_rpm*2*pi/60; % [rad s�?1]
 y0 = [sub_theta2, sub_dtheta2];


 % Time Integration
 numbofrev = 2;
 tint = abs(numbofrev/(sub_dtheta2_rpm/60));
 dtint = 1e-3;


 tol = 1e-12;
maxit = 10;

% 1
sub_init = CalcAngles(y0, lengths);


syms m g l gamma M I
syms theta2 theta3 theta4 theta5
syms dtheta2 dtheta3 dtheta4 dtheta5
syms ddtheta2 ddtheta4 ddtheta5
syms O2A O4B BC O4O2 O4G4 BG5 yc m3 m4 m5 m6 J4 J5 J2
syms lambda1 lambda2

%% values
O2A = 0.2;
O4B = 0.7;
BC = 0.6;
O4O2 = 0.3;
O4G4 = 0.4;
BG5 = 0.3;
yC = 0.9;
m3 = 0.5;
m4 = 6;
m5 = 4;
m6 = 2;
J2 = 100;
J4 = 10;
J5 = 6;
theta3 = theta4;

%% Gen coord

q = [theta2; theta4; theta5];
qd = [dtheta2; dtheta4; dtheta5];
qdd = [ddtheta2; ddtheta4; ddtheta5];

Q = [0;0];

%% Locations

x3 = O2A*cos(theta2);
y3 = O4O2 + O2A*sin(theta2);

x4 = O4G4*cos(theta4);
y4 = O4G4*sin(theta4);

x5 = O4B*cos(theta4) + BG5*cos(theta5);
y5 = O4B*sin(theta4) + BG5*sin(theta5);

x6 = O4B*cos(theta4) + BC*cos(theta5);
y6 = O4B*sin(theta4) + BC*sin(theta5);

%% Velocities

Ti = [theta2; x3;y3;x4;y4;theta4;x5;y5;theta5;x6];
Ti_k = jacobian(Ti,q);

% velocity
dTi = jacobian(Ti,q)*qd;

ddTi = jacobian(Ti,q)*qdd;
K = [J2 m3 m3 m4 m4 J4 m5 m5 J5 m6];
Mij = diag(K);



fi = [0;0;0;0;0;0;0;0;0;1000];
hj = jacobian(Ti_k*qd,q)*qd;

%% Constraints

l = sqrt(x3^2 + y3^2);
Di = [O4B*sin(theta4) + BC*sin(theta5) - yc;O2A*cos(theta2) - l*cos(theta4)];
lambda = [lambda1;lambda2];
Di_j = jacobian(Di,q);
Di_jk = jacobian((Di_j*qd),q);
hj_constraint = -jacobian(Di_j*qd,q)*qd;

%% TMT solution

T = Ti_k;
M = Mij;
G = hj;
dT = dTi;
ddT = ddTi;
Mred = T.'*M*T;
Mred = simplify(Mred);
pretty(Mred);
fred = T.'*fi - T.'*M*G;
fred = simplify(fred);
pretty(fred);

A = [Mred Di_j.'; Di_j zeros(2,2)];
B = [fred; hj_constraint];

sol = A\B;

 
% 2 

y_fin = Integrator(sub_init, tint, dtint, tol, maxit);


[row_y_fin, col_y_fin] = size(y_fin);
time_stepsize = (row_y_fin*dtint)/(row_y_fin - 1);
t_new = 0:time_stepsize:(row_y_fin*dtint);
if exist('vel3', 'var') == 0
    vel3 = zeros(row_y_fin, 1);

    for n = 1:row_y_fin
 theta2 = y_fin(n, 1);
 dtheta2 = y_fin(n, 4);
 vel3x = eval(dT(2));
 vel3y = eval(dT(3));
 vel3(n,1) = cos(theta4)*vel3x + sin(theta4)*vel3y;
    end
else
end

% slider 6 + constraint forces

if exist('states6', 'var') == 0
states6 = zeros(row_y_fin, 3);
multipliers = zeros(row_y_fin, 2);

 for n = 1:row_y_fin
 theta2 = y_fin(n, 1);
 dtheta2 = y_fin(n, 4);
 theta4 = y_fin(n, 2);
 dtheta4 = y_fin(n, 5);
 theta5 = y_fin(n, 3);
 dtheta5 = y_fin(n, 6);

 [dqandddq, Multipliers] = DAE(0, [theta2, theta4, theta5, dtheta2, dtheta4, dtheta5]);
 ddtheta4 = dqandddq(5);
 ddtheta5 = dqandddq(6);
 states6(n, :) = [eval(T(10)), eval(dT(10)), eval(ddT(10))];
 multipliers(n, :) = Multipliers;
 end
else
    end

%% PLOTS
figure
 plot(t_new, y_fin(:, 4), 'r-');
 hold on
 plot(t_new,y_fin(:, 5), 'b-');
 plot(t_new, y_fin(:, 6), 'k-');
 xlabel('Time [s]')
 ylabel('Angular velocity [rads/s]')
 title('Angular Velocities')
 inf_question_b = legend('$Crank({\theta}_{2})$', '$Rocker({\theta}_{4})$','$Connecting Bar({\theta}_{5}$)');

 set(inf_question_b, 'Interpreter', 'Latex', 'FontSize', 14);

hold off

figure
 plot(t_new, vel3)
 xlabel('Time [sec]')
 ylabel('Speed m/s')
 title('Velocity of slider 3 wrt rocker 4')
 inf_question_c = legend('Relative Velocity')

 figure
 plot(t_new, states6(:, 1))
 xlabel('Time [s]')
 ylabel('Position of Slider 6 [m]')
 title('Position of slider 6')
 figure
 plot(t_new, states6(:, 2))
 xlabel('Time [s]')
 ylabel('Velocity of slider 6 [m/s]')
 title('Velocity of slider 6')
 figure
 plot(t_new, states6(:, 3))
 xlabel('Time [s]')
 ylabel('Acceleration of slider 6 [m/s]')
 title('Acceleration of slider 6')
 
figure
 plot(t_new, multipliers(:,1), '--')
 hold on
plot(t_new, multipliers(:,2), '-.')
 xlabel('Time [s]')
 ylabel('Normal Force [N]')
 title('Normal forces of slider 6 on ground and slider 3 on rocker 4')
 inf_question_d = legend(' Slider 3 on rocker 4', 'Slider 6 on the ground');
 set(inf_question_d, 'FontSize', 14);
hold off

%% DAE 

function [SolutionAngles, Multipliers] = DAE(t,y)
theta2 = y(1);
theta4 = y(2);
theta5 = y(3);
dtheta2 = y(4);
dtheta4 = y(5);
dtheta5 = y(6);


A = [5001/50,0,0,0,(cos(theta4)*((2*cos(theta2)*sin(theta2))/25 - (2*cos(theta2)*(sin(theta2)/5 + 3/10))/5))/(2*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2)) - sin(theta2)/5;
0,(49*sin(theta4)^2)/50 + 323/25,(63*cos(theta4 - theta5))/50 - (21*cos(theta4 + theta5))/50, (7*cos(theta4))/10,sin(theta4)*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2);
0,(63*cos(theta4 - theta5))/50 - (21*cos(theta4 + theta5))/50,  (18*sin(theta5)^2)/25 + 159/25,  (3*cos(theta5))/5,                                                                                                                                                                 0;
                                                                                                                                                                 0,                                              (7*cos(theta4))/10,                                           (3*cos(theta5))/5,                  0,                                                                                                                                                                 0;
(cos(theta4)*((2*cos(theta2)*sin(theta2))/25 - (2*cos(theta2)*(sin(theta2)/5 + 3/10))/5))/(2*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2)) - sin(theta2)/5, sin(theta4)*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2),                                                           0,                  0,                                                                                                                                                                 0];



B = [0;...
 (21*dtheta5^2*cos(theta4)*sin(theta5))/25 - (49*dtheta4^2*cos(theta4)*sin(theta4))/50 - 700*sin(theta4) - (42*dtheta5^2*cos(theta5)*sin(theta4))/25;...
 (21*dtheta4^2*cos(theta5)*sin(theta4))/25 - (42*dtheta4^2*cos(theta4)*sin(theta5))/25 - 600*sin(theta5) - (18*dtheta5^2*cos(theta5)*sin(theta5))/25;...
 (7*sin(theta4)*dtheta4^2)/10 + (3*sin(theta5)*dtheta5^2)/5;...
 dtheta2*(dtheta2*(cos(theta2)/5 - (cos(theta4)*((2*cos(theta2)*sin(theta2))/25 - (2*cos(theta2)*(sin(theta2)/5 + 3/10))/5)^2)/(4*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(3/2)) + (cos(theta4)*((2*sin(theta2)^2)/25 - (2*sin(theta2)*(sin(theta2)/5 + 3/10))/5))/(2*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2))) + (dtheta4*sin(theta4)*((2*cos(theta2)*sin(theta2))/25 - (2*cos(theta2)*(sin(theta2)/5 + 3/10))/5))/(2*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2))) - dtheta4*(dtheta4*cos(theta4)*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2) - (dtheta2*sin(theta4)*((2*cos(theta2)*sin(theta2))/25 - (2*cos(theta2)*(sin(theta2)/5 + 3/10))/5))/(2*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2)))];
ddqandmultipliers = A\B;
SolutionAngles = [dtheta2; dtheta4; dtheta5; ddqandmultipliers(1:3)];
Multipliers = ddqandmultipliers(4:5);
end
%% Angle calculation

function sub_init = CalcAngles(y, lengths)

syms init_theta2 real
syms dtheta2_0

O2A = lengths(1);
O4B = lengths(2);
BC = lengths(3);
O4O2 = lengths(4);
yC = lengths(7);

% Initial pos

init_theta4 = atan((O4O2 + O2A*sin(init_theta2))/(O2A*cos(init_theta2)));
theta5_0 = (pi/2)+acos((yC - O4B*sin(init_theta4))/BC);

% init vel
dtheta4_0 = jacobian(init_theta4,init_theta2)*dtheta2_0;
dtheta5_0 = jacobian(theta5_0,init_theta2)*dtheta2_0;

% init condition

theta2_s = y(1);
dtheta2_s = y(2);

% Substitution

theta4_s = subs(init_theta4,init_theta2, theta2_s);
theta5_s = subs(theta5_0, init_theta2, theta2_s);
dtheta4_s = subs( dtheta4_0, [init_theta2, dtheta2_0], [theta2_s dtheta2_s]);
dtheta5_s = subs(dtheta5_0, [init_theta2, dtheta2_0], [theta2_s dtheta2_s]);

sub_init = eval([theta2_s; theta4_s; theta5_s; dtheta2_s; dtheta4_s; dtheta5_s]);

end
%% COnstraints

function [ D, dD] = Constraints(phi)

theta2 = phi(1);
theta4 = phi(2);
theta5 = phi(3);

O2A = 0.2;
O4B = 0.7;
BC = 0.6;
O4O2 = 0.3;
O4G4 = 0.4;
BG5 = 0.3;
yc = 0.9;
m3 = 0.5;
m4 = 6;
m5 = 4;
m6 = 2;
J2 = 100;
J4 = 10;
J5 = 6;
theta3 = theta4;

x3 = O2A*cos(theta2);
y3 = O4O2 + O2A*sin(theta2);

x4 = O4G4*cos(theta4);
y4 = O4G4*sin(theta4);

x5 = O4B*cos(theta4) + BG5*cos(theta5);
y5 = O4B*sin(theta4) + BG5*sin(theta5);

x6 = O4B*cos(theta4) + BC*cos(theta5);
y6 = O4B*sin(theta4) + BC*sin(theta5);
l = sqrt(x3^2 + y3^2);
D = [O4B*sin(theta4) + BC*sin(theta5) - yc;O2A*cos(theta2) - l*cos(theta4)];
dD = [ 0,(7*cos(theta4))/10 , (3*cos(theta5))/5 ;...
    (cos(theta4)*((2*cos(theta2)*sin(theta2))/25 - (2*cos(theta2)*(sin(theta2)/5 + 3/10))/5))/(2*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2)) - sin(theta2)/5, sin(theta4)*((sin(theta2)/5 + 3/10)^2 + cos(theta2)^2/25)^(1/2),0];
end

%% Integration


function y_fin = Integrator(y,t,dt,tol,maxit)
step = 0;
theta2 = y(1);
y_fin(1,:) = y.';

% RK4 Method

while abs(theta2) < (4*pi)
    step = step+1;
    [k1,~] = DAE(0,y);
    y2 = y + (dt/2).*k1;
    [k2,~] = DAE(0,y2);
    y3 = y + (dt/2).*k2;
    [k3,~] = DAE(0,y3);
    y4 = y + (dt).*k3;
    [k4,~] = DAE(0,y4);
    y = y+ (dt/6).*(k1+2*k2+2*k3+k4);
    
    y_t = y(1:3);
    [D,dD] = Constraints(y_t);
    numbit = 0;
    
    while ((max(abs(D))> tol) && (numbit < maxit))
        
        numbit = numbit +1;
        delta_q = -dD.'*inv(dD*dD.')*D;
        y_t = y_t + delta_q;
        [D,dD] = Constraints(y_t);
    end
    
    dy_t = y(4:6);
    [~,dD] = Constraints(y_t);
    delta_dq = -dD.'*inv(dD*dD.')*dD*dy_t;
    dy_t = dy_t + delta_dq;
    
    y = [y_t; dy_t];
    theta2 = y_t(1,1);
    y_fin(step+1,:) = y.';
end
end
 