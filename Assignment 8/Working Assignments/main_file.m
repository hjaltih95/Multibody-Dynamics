clc
clear
close all

global A_ M1_ M2_ I1_ I2_ F_ FC_ B_ C_ D_

%% set globals and initial conditions
A_ = 0.2; % distance skate A to cm 1
B_ = 0.6; % distance cm 1 to hinge B
C_ = 0.1; % distance cm 2 to skate C
D_ = 0.1; % distance cm 2 to hinge B
M1_ = 46; % mass sleigh 1
M2_ = 1; % mass sleigh 2
I1_ = 10; % mass moment of inertia of sleigh at cm 1
I2_ = 0.005; % mass moment of inertia of sleigh at cm 2

M = [M1_ M1_ I1_ M2_ M2_ I2_];

% constant applied forces

F_ = [0,0,0,0,0,0].';

%Initial conditions for first part uncomment when needed

% x01 = A_; % [m]
% y01 = 0; % [m]
% phi01 = 0; % [rad]
% phi02 = pi/2; % [rad]
% x02 = A_ + B_; % [m]
% y02 = D_; % [m]
% xd01 = 1; % [m/s]
% yd01 = 0; % [m/s]
% phid01 = 0; % [rad/s]
% xd02 = 0.5; % [m]
% yd02 = 0; % [m]
% phid02 = 5; % [rad/s]

% Initial Conditions for 2nd part

x01 = 0; % [m]
y01 = 0; % [m]
phi01 = 0; % [rad]
phi02 = pi; % [rad]
x02 = 0 ; % [m]
y02 = 0; % [m]
xd01 = 0; % [m/s]
yd01 = 0; % [m/s]
phid01 = 0; % [rad/s]
xd02 = 0; % [m]
yd02 = 0; % [m]
phid02 = 0; % [rad/s]

% initial conditions in a vector
y0 = [x01 y01 phi01 x02 y02 phi02 xd01 yd01 phid01 xd02 yd02 phid02].';
t0 = 0;
tend = 60;
% stepsize
dt = 0.01;
%% simulate motion
nstep = (tend-t0)/dt+1;
t = t0;
y = y0;
T = zeros(nstep,1);
Y = zeros(nstep,12);
Ev = zeros(nstep,1);
Es = zeros(nstep,1);
T(1) = t;
Y(1,:) = y.';
for i = 2:nstep
f = AppliedForces(t,y)
y_start = y;
[t,y]=rk4step(@odedaeEzy,t,y,dt);
% project the positions on the constraint surface
y = projectPosition(t,y);
% project the speeds on the constraint surface
y = projectSpeed(t,y);
T(i) = t;
Y(i,:) = y.';

%Torque
M1(i) = 2*cos(pi*T(i));
end


%% post processing
X1 = Y(:,1) ; Y1 = Y(:,2) ; P1 = Y(:,3);
X2 = Y(:,4); Y2 = Y(:,5); P2 = Y(:,6);
DX1 = Y(:,7) ; DY1 = Y(:,8) ; DP1 = Y(:,9);
DX2 = Y(:,10); DY2 = Y(:,11); DP2 = Y(:,12);
XA = X1- A_*cos(P1) ; YA = Y1 - A_*sin(P1);
XB1 = X1+ B_*cos(P1) ; YB1 = Y1+B_*sin(P1);
XB2 = X2 - D_*cos(P2) ; YB2 = Y2-D_*sin(P2);
XC = X2+ C_*cos(P2) ; YC = Y2 + C_*sin(P2);
% Linear velocity for body 1
LV1 = sqrt(DX1.^2 + DY1.^2);
% Relative angular velocity
RVB = DP2 - DP1;

% Kinetic Energy
K_kin1 = 0.5* M1_ * (DX1.^2 + DY1.^2) + 0.5* I1_ *DP1.^2;
K_kin2 = 0.5* M2_ * (DX2.^2 + DY2.^2) + 0.5* I2_ *DP2.^2;

E_ke = K_kin1 + K_kin2;

% Work done
for j = 1 :length(T)
    if j ==1
        W(j) = 0;
    else
        W(j) = (2*cos(pi*T(j)))*(P1(j) - P1(j-1)) - (2*cos(pi*T(j)))*(P2(j) - P2(j-1));
    end
end

% E Path of COM of BODY 1 and BODY 2
figure
plot(X1,Y1)
hold on
grid on
plot(X2,Y2)
xlabel('X Position[m]')
ylabel('Y position [m]')
title('COM of Body 1 and Body 2')
legend('COM Body 1','COM Body 2')

%F) Plot the path of A and C
figure
plot(XA,YA,XC,YC)
grid on
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point A','Point C')
title('Path of A and C')

% G LINEAR AND ANGULAR VEL OF BODY 1 WRT TIME
figure
grid on
subplot(3,1,1)
plot(T,LV1)
title('Linear Velocity')
subplot(3,1,2)
plot(T,DX1)
title('Angular Velocity')
ylabel('Linear velocities [m/s] and angular velocies [rad/s]')
subplot(3,1,3)
plot(T,LV1)
hold on
plot(T,DX1)
xlabel('Time [s]')
legend('Linear Velocity','Angular Velocity');

% H) KE of system
figure
plot(T,E_ke)
xlabel('Time [s]')
ylabel('Energy [J]')
legend('Kinetic energy of the system');
title('Kinetic Energy of the system wrt Time')

% I Torque and rel angular velocity
figure
plot(T,M1)
hold on
grid on
plot(T, RVB)
xlabel('Time [s]');
ylabel('Torque[N-m] and Rel Angular Vel[rad/s]');
legend('Torque','Rel Vel');
title('Torque and Relative Angular Velocity');

%J) Kinetic energy of the system and work done by the torque
figure
plot(T,E_ke)
hold on
grid on
plot(T, cumsum(W))
xlabel('Time [s]')
ylabel('Energy [J]')
legend('Kinetic energy of the system','Work')
