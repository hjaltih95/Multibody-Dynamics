clc;clear;

%% setup

syms x1 y1 x2 y2 phi1 phi2;
syms dx1 dy1 dx2 dy2 dphi1 dphi2;
syms ddx1 ddy1 ddx2 ddy2 ddphi1 ddphi2;
%syms a b c d;
%syms m1 m2 I1 I2

global a b c d m1 m2 I1 I2 

a = 0.2; b = 0.6; c = 0.1; d = 0.1; m1 = 46; m2 = 1; I1 = 10; I2 = 0.005;


x = [x1; y1; phi1; x2; y2; phi2];
dx = [dx1; dy1; dphi1; dx2; dy2; dphi2];

% holonomic constraints
C1 = x1 + b*cos(phi1) - x2 + d*cos(phi2);
C2 = y1 + b*sin(phi1) - y2 + d*sin(phi2);
C = [C1; C2];

Cd = jacobian(C, x);
Cd = simplify(Cd);
% convective terms
Ch = jacobian(Cd*dx, x)*dx;
Ch = simplify(Ch);

% nonholonomic constraints
S1 = -dx1*sin(phi1) + dy1*cos(phi2) - a*dphi1;
S2 = -dx2*sin(phi2) + dy2*cos(phi2) + c*dphi2;
S = [S1; S2];

% convective terms
Sh = [-dphi1*dx1*cos(phi1)-dphi1*dy1*sin(phi1);
      -dphi2*dx2*cos(phi2)-dphi2*dy2*sin(phi2)];
  
% DAE
M = diag([m1 m1 I1 m2 m2 I2]);
Sd = [-sin(phi1) 0; cos(phi1) 0; -a 0; 0 -sin(phi2); 0 cos(phi2); 0 c];
A = [M Cd.' Sd; Cd zeros(2,4); Sd.' zeros(2,4)];
B = [zeros(6,1); -Ch; -Sh];
Acc = simplify(A\B);
matlabFunction(Acc,'File','Acceleration');
matlabFunction(Cd, 'File', 'HoC');
matlabFunction(Ch, 'File', 'HoCc');
matlabFunction(Sd, 'File', 'nHoC');
matlabFunction(Sh, 'File', 'nHoCc');



%% Initial conditions
% x1 = a;
% y1 = 0;
% phi1 = 0; 
% x2 = a+b;
% y2 = d;
% phi2 = pi/2;
% dx1 = 1; 
% dphi1 = 0;
% 
% 
% S1 = -dx1*sin(phi1) + dy1*cos(phi2) - a*dphi1;
% S2 = -dx2*sin(phi2) + dy2*cos(phi2) + c*dphi2;
% Cd = HoC(phi1,phi2);
% dx = [dx1; dy1; dphi1; dx2; dy2; dphi2];
% 
% Cons = Cd*dx;
% 
% equ1 = Cons(1);
% equ2 = Cons(2);
% equ3 = S1;
% equ4 = S2;
% sol = solve(equ1==0, equ2==0, equ3==0, equ4==0);
% 
% % solved initial conditions
% dx2 = vpa(sol.dx2);
% dy2 = round(vpa(sol.dy2));
% dy1 = vpa(sol.dy1);
% dphi2 = vpa(sol.dphi2);
% 
% % initial conditions
% y0 = [x1, y1, phi1, x2, y2, phi2, dx1, dy1, dphi1, dx2, dy2, dphi2]';
% 
% % time
% t0 = 0;
% t_end = 1;
% dt = 0.01;
%% Simulation
% numsteps = (t_end - t0)/dt + 1;
% t = t0;
% y = y0;
% % create empty vectors + first values
% T = zeros(numsteps, 1);
% Y = zeros(numsteps, 12);
% T(1) = t;
% Y(1,:) = y.';
% 
% for i = 2:numsteps
% % Numerical integration by means of RK 4th order
% [t,y] = rk4step(@odedae,t,y,dt);
% % Coordinate projection Method
% y = projectSpeed(t,y);
% T(i) = t;
% Y(i,:) = y.';
% end
% 
% % Declaration of values
% X1 = Y(:,1); Y1 = Y(:,2); PHI1 = Y(:,3);
% X2 = Y(:,4); Y2 = Y(:,5); PHI2 = Y(:,6);
% X1D = Y(:,7); Y1D = Y(:,8); PHI1D = Y(:,9);
% X2D = Y(:,10); Y2D = Y(:,11); PHI2D = Y(:,12);
% XA = X1 - a * cos(PHI1);
% YA = Y1 - a * sin(PHI1);
% XC = X2 + c * cos(PHI2);
% YC = Y2 + c * sin(PHI2);
% 
% % Plot
% % Path question D
% figure
% plot(X1,Y1)
% hold on
% %plot(XA,YA)
% %hold on
% plot(X2,Y2)
% hold on
% %plot(XC,YC)
% %axis equal
% title('QuestionD - Path of the CoM')
% legend('CM1','CM2')
% xlabel('x [m]')
% ylabel('y [m]')
% grid on


%% Question E
% Initial conditions

x1 = 0;
y1 = 0;
phi1 = 0;
phi2 = pi;
dx1 = 0;
dy1 = 0;
dphi1 = 0;
dphi2 = 0;
x2 = x1 + b * cos(phi1) + d * cos(phi2);
y2 = y1 + b * sin(phi1) + d * sin(phi2);
dx2 = dx1 - b * dphi1 * sin(phi1) - d * dphi2 * sin(phi2);
dy2 = dy1 + b * dphi1 * cos(phi1) + d * dphi2 * cos(phi2);
y0 = [x1, y1, phi1, x2, y2, phi2, dx1, dy1, dphi1, dx2, dy2, dphi2].';
t0 = 0;
tend = 60;
td = 0.005; % Stepsize

% Additional values
M0 = 2;
omega = pi;
M1 = M0 * cos(omega *td);
F = [0, 0, -M1, 0, 0, M1].';
nstep = (tend - t0) / td +1 ;
t = t0;
y = y0;
T = zeros(nstep,1);
Y = zeros(nstep,12);
W = zeros(nstep,1);
Es = zeros(nstep,1);
T(1) = t;
Y(1,:) = y.';




for i = 2 : nstep
[t,y] = rk4step(@odedae2,t,y,td);
 y_start = y0;
 y = projectSpeed(t,y);
 %y = projectposition(y);
f = AppliedForcesE(t,y);
T(i) = t;
Y(i,:) = y.';
end

% Post Processing
X1 = Y(:,1); Y1 = Y(:,2); PHI1 = Y(:,3);
X2 = Y(:,4); Y2 = Y(:,5); PHI2 = Y(:,6);
X1D = Y(:,7); Y1D = Y(:,8); PHI1D = Y(:,9);
X2D = Y(:,10); Y2D = Y(:,11); PHI2D = Y(:,12);
XA = X1 - a * cos(PHI1);
YA = Y1 - a * sin(PHI1);
XC = X2 + c * cos(PHI2);
YC = Y2 + c * sin(PHI2);
for j = 1:length(T)
    if j==1
        W(j) = 0;
    else
        W(j) = 2*cos(pi*T(j))*(PHI1(j)-PHI1(j-1))-2*cos(pi*T(j))*(PHI2(j)-PHI2(j-1));
    end
end

% Kinetic energy
E_kin1 = 1/2 * m1 * (X1D.^2 + Y1D.^2) + 1/2 * I1 * PHI1D.^2;
E_kin2 = 1/2 * m2 * (X2D.^2 + Y2D.^2) + 1/2 * I2 * PHI2D.^2;
E_kin = E_kin1 + E_kin2;
% Plot
figure(1)
plot(X1,Y1)
hold on
% plot(XA,YA)
% hold on
plot(X2,Y2)
% hold on
% plot(XC,YC)
axis equal
title('Question f - Path of Skates A and C together with the CM')
legend('A','C')
xlabel('x [m]')
ylabel('y [m]')
grid on

figure(2)
speed = sqrt(X1D.^2 + Y1D.^2);
plot(T,speed)
hold on
plot(T,PHI1D)
title('Qg - Speed and velocity of body 1')

figure(3)
M=2*cos(omega*T);
plot(T,M)
hold on
plot(T,PHI2D - PHI1D)
xlabel('t [sec]')
ylabel('X\phiD/TD')
legend('\phi1d/td','\phi2d/td')
title('Qi - Torque and relative velocity')
grid on

figure(4)
plot(T, E_kin)
hold on
plot(T, cumsum(W))
title('Question h and j - Work and Kinetic Energy')
xlabel('t [sec]')
ylabel('Energy [J]')
legend('EKinetic','Work')

%%
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
    
    delta_p = -Dv.'*((Dv*Dv.')\s);
    q_corr(1:6) = q_corr(1:6) + delta_p;
    
    i = i+1;
    
end

end