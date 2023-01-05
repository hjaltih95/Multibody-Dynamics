close all
clear all
clc
delete hw7_x6Function.m hw7_constraintFunction.m hw7_qddotFunction.m

%% parameters
%lenghts
l2 = 0.2;
l4 = 0.7;
l5 = 0.6;
l4G = 0.4;                                                                  %distance to center of gravity
l5G = 0.3;                                                                  %distance to center of gravity
yC = 0.9;                                                                   %height of the slider
yO2 = 0.3;                                                                  %height of the crank center

%masses
m3 = 0.5;
m4 = 6;
m5 = 4;
m6 = 2;

%moments of inertia
J4 = 10;
J5 = 6;
J2 = 100;

%applied forces and torques
F = 1000;
T = 0;

%% equations of motion
syms theta2 theta4 theta5 theta2dot theta4dot theta5dot theta2ddot theta4ddot theta5ddot t real

%generalized coordinates, velocities and accelerations
q = [theta2 theta4 theta5].';
qdot = [theta2dot theta4dot theta5dot].';
qddot = [theta2ddot theta4ddot theta5ddot].';                               %accelerations are only required for the calculation of x6dot
yhalf = length(q);
yhalfplus1 = length(q)+1;

%CoM coordinates (velocities and accelerations not required for TMT)
x = [l2*cos(theta2) yO2+l2*sin(theta2) theta2...                            %x = [x3 y3 theta2...
    l4G*cos(theta4) l4G*sin(theta4) theta4...                               %x4 y4 theta4...
    l4*cos(theta4)+l5G*cos(theta5) l4*sin(theta4)+l5G*sin(theta5) theta5... %x5 y5 theta5...    
    l4*cos(theta4)+l5*cos(theta5) l4*sin(theta4)+l5*sin(theta5)].';         %x6 y6].' 

%x-coordinate, x-velocity and x-acceleration of slider 6
x6 = x(10);
x6dot = simplify(jacobian(x6,q)*qdot);
x6ddotval = simplify(jacobian(x6dot,q)*qdot + jacobian(x6dot,qdot)*qddot);
matlabFunction([x6, x6dot, x6ddotval],'file','hw7_x6Function','vars',{[q; qdot; qddot]}); %save the equations to a function

%generalized mass matrix
M = diag([m3,m3,J2,m4,m4,J4,m5,m5,J5,m6,m6]);
Tg = simplify(jacobian(x,q));                                               %transformation matrix from CoM to generalized coordinates
Mg = simplify(Tg.'*M*Tg);

%additional slider constraints for the calculation of the normal forces
C3 = simplify(x(2)/x(1)-x(5)/x(4));                                         %constraint for slider 3
C6 = simplify(yC - x(11));                                                  %constraint for slider 6
C = [C3; C6];                                                               %constraint matrix
Cq = simplify(jacobian(C,q));                                               %c,q
Cqqqdot2 = simplify(jacobian(Cq*qdot,q)*qdot);                              %c,qq*qdot*qdot
matlabFunction(C,Cq,'file','hw7_constraintFunction','vars',{q});            %save the equations to a function

%applied forces on generalized coordinates (generalized forces)
Q = zeros(length(q),1);

%applied forces on CoM coordinates
F = [0 0 0 0 0 0 0 0 0 F 0].';

%transformation of velocities and accelerations
qdotT = Tg*qdot;                                                            %transformed velocities
qddotT = simplify(jacobian(qdotT,q)*qdot);                                  %transformed accelerations

%left- and right-hand side
lhs = simplify([Mg Cq.'; Cq zeros(length(C))]);
rhs = simplify([Q + Tg.'*(F-M*qddotT); -Cqqqdot2]);
qddot = simplify(lhs\rhs);
matlabFunction(qddot,'file','hw7_qddotFunction','vars',{[q; qdot]});        %save the equations to a function

%% calculations of initial conditions from theta2 and theta2dot
%provided conditions theta2, theta2dot
theta2val0 = 0;                                                             
theta2dotval0 = 75/60*2*pi;

%theta4, theta4dot
x3sym0 = x(1);               
y3sym0 = x(2);
theta4sym0 = atan(y3sym0/x3sym0);                                           %tan(theta4) = y3(phi2)/x3(phi2)
theta4dotsym0 = jacobian(theta4sym0,theta2)*theta2dot;

%theta5, theta5dot
yBC = yC-l4*sin(theta4sym0);
theta5sym0 = pi-asin(yBC/l5);                                               %sin(pi-theta5) = yBC/l5
theta5dotsym0 = jacobian(theta5sym0,theta2)*theta2dot;

%complete initial condition vector
ysym0 = [theta2 theta4sym0 theta5sym0 theta2dot theta4dotsym0 theta5dotsym0].';     %how to calculate the initial conditions
yval0 = double(subs(ysym0, [theta2; theta2dot], [theta2val0; theta2dotval0]));      %calculate the initial conditions

%% simulation

tol = 10^-12;                                                               %tolerance for coordinate projection
maxiterat = 10;                                                             %maximum iterations for coordinate projection
count = 0;     
te = 4*pi/theta2dotval0+0.2;                                                %simulation duration (a little over two revolutions)
n = 14;                                                                     %2^n steps
h = te/(2^n);                                                               %step size in seconds
N = ceil(te/h);                                                             %number of steps until te
yn0 = yval0;   

for i = 1:N                                                                 %number of steps
    
    %RK4-integration
    ydot0 = hw7_qddotFunction(yn0);
    k1 = [yn0(yhalfplus1:end);ydot0(1:yhalf)];
    yn1 = yn0 + 1/2*h*k1;
    ydot1 = hw7_qddotFunction(yn1);
    k2 = [yn1(yhalfplus1:end);ydot1(1:yhalf)];
    yn2 = yn0 + 1/2*h*k2;
    ydot2 = hw7_qddotFunction(yn2);
    k3 = [yn2(yhalfplus1:end);ydot2(1:yhalf)];
    yn3 = yn0 + h*k3;
    ydot3 = hw7_qddotFunction(yn3);
    k4 = [yn3(yhalfplus1:end);ydot3(1:yhalf)];
    yn4 = yn0 + 1/6*h*(k1+2*k2+2*k3+k4);
    yn1 = yn4;                                                              %y_n+1 before coordinate projection    
    
    %Gauss-Newton coordinate projection (multiple iterations)
    iterat = 0; 
    [Cn, Cnq] = hw7_constraintFunction(yn1(1:yhalf));                       %set iteration for coordinate projection to 0
    while max(max(abs(Cn))) > tol && iterat <= maxiterat
        Dqn1 = -Cnq.'*inv(Cnq*Cnq.')*Cn;                                    %smallest (orthogonal) distance between calculated q_n+1 and the constraint surface
       	yn1(1:yhalf) = yn1(1:yhalf) + Dqn1;                                 %q_n+1 after coordinate projection
        [Cn, Cnq] = hw7_constraintFunction(yn1(1:yhalf));
        iterat = iterat+1;
    end      
    
    %Gauss-Newton velocity projection (single iteration, linear problem)
    Dqdotn1 = -Cnq.'*inv(Cnq*Cnq.')*Cnq*yn1(yhalfplus1:end);
   	yn1(yhalfplus1:end) = yn1(yhalfplus1:end) + Dqdotn1;                    %qdot_n+1 after coordinate projection
    
    %closing of iteration
    yval(i,:) = yn1;                                                        %y for each time step
    ydotval(i,:) = hw7_qddotFunction(yn1);                                  %ydot for each time step
    tval(i) = i*h;                                                          %timestamps for each step
    yn0 = yn1;                                                              %set the next coordinates and velocities to be the current coordinates and velocities
    count = count+1;                                                        %add 1 to the iteration count
    disp(num2str(round(count/N*100,2)+"% - RK4"))                           %display the progress percentage
end

%derivation of the required parameters from q_n+1 and qdot_n+1
qval = yval(:,1:yhalf);
qdotval = yval(:,yhalfplus1:end);
qddotval = ydotval(:,1:yhalf);
lagval = ydotval(:,yhalfplus1:end);                                         %Lagrangian multipliers
v34 = l2.*qdotval(:,1).*cos(qval(:,1)-qval(:,2)+pi/2);                      %velocity of slider 3 with respect to bar 4
F34 = lagval(:,1);                                                          %constraint force of slider 3, perpendicular to bar 4
F6 = lagval(:,2);                                                           %vertical constraint force of slider 6
for i = 1:length(tval)
    x6vals(i,:) = hw7_x6Function([qval(i,:) qdotval(i,:) qddotval(i,:)]');
end
x6val = x6vals(:,1);
x6dotval = x6vals(:,2);
x6ddotval = x6vals(:,3);

%calculate period
[~, Tp] = findpeaks(yval(:,5));
period = (Tp(2)-Tp(1))*h;

%% plot results
close all

%theta4dot, theta5dot, theta6dot
figure; hold on
plot(tval,qdotval(:,1))
plot(tval,qdotval(:,2))
plot(tval,qdotval(:,3))
xlabel('time [s]')
ylabel('angular velocity [rad/s]')
title('angular velocity of crank 2, rocker 4 and connecting bar 5')
legend('\theta_2','\theta_4','\theta_5','Location','southwest')

%velocity of slider 3 with respect to bar 4
figure
plot(tval,v34)
xlabel('time [s]')
ylabel('velocity [m/s]')
title('sliding velocity of slider 3 with respect to rocker 4')
legend('v_{3-4}')

%normal forces F34 and F6
figure; hold on
plot(tval,F6)
plot(tval,F34)
xlabel('time [s]')
ylabel('force [N]')
title('normal forces of sliders 3 and 6')
legend('F_3','F_6','Location','southwest')

%x-coordinate of slider 6
figure
plot(tval,x6val)
xlabel('time [s]')
ylabel('position [m]')
title('x-coordinate of slider 6.')
legend('x_6','Location','northwest')

%x-velocity of slider 6
figure
plot(tval,x6dotval)
xlabel('time [s]')
ylabel('velocity [m/s]')
title('x-velocity of slider 6.')
legend('v_6','Location','northwest')

%x-acceleration of slider 6
figure
plot(tval,x6ddotval)
xlabel('time [s]')
ylabel('acceleration [m/s^2]')
title('x-acceleration of slider 6.')
legend('a_6','Location','northwest')

%% animation
close all
%calculate joint coordinates
xAi = l2.*cos(qval(:,1));
yAi = yO2+l2.*sin(qval(:,1));
xG4i = l4G.*cos(qval(:,2));
yG4i = l4G.*sin(qval(:,2));
xBi = l4.*cos(qval(:,2));
yBi = l4.*sin(qval(:,2));
xG5i = l4.*cos(qval(:,2))+l5G.*cos(qval(:,3));
yG5i = l4.*sin(qval(:,2))+l5G.*sin(qval(:,3));
xCi = l4.*cos(qval(:,2))+l5.*cos(qval(:,3));
yCi = l4.*sin(qval(:,2))+l5.*sin(qval(:,3));

step = 30;                                                                  %number of frames to skip between steps
tstop = 1;                                                                  %animation end-time (simulation time)
if tstop <= te
    indexstop = round(tstop/te*length(tval));
else
    indexstop = length(tval);
end

%skip steps to increase speed
xAi = xAi(1:step:indexstop);
yAi = yAi(1:step:indexstop);
xG4i = xG4i(1:step:indexstop);
yG4i = yG4i(1:step:indexstop);
xBi = xBi(1:step:indexstop);
yBi = yBi(1:step:indexstop);
xG5i = xG5i(1:step:indexstop);
yG5i = yG5i(1:step:indexstop);
xCi = xCi(1:step:indexstop);
yCi = yCi(1:step:indexstop);

%animate
xlimval = [-1 0.5];
ylimval = [0 1];
fixsize = 5;
rsize = 5;
ssize = 8;
for i = 1:length(xAi)
    clf
    hold on; axis equal
    
    %plot rotary joints
    pB = plot(xBi(i),yBi(i),'or','MarkerSize',rsize)
    
    %plot slider joints and slider lines
    p3 = plot(xAi(i),yAi(i),'sr','MarkerSize',ssize);
    p6 = plot(xCi(i),yCi(i),'sr','MarkerSize',ssize);
    p6l = plot([xlimval(1),xlimval(2)],[yC,yC],':k');
    
    %plot fixed joints
    pO2 = plot(0,yO2,'ok','MarkerSize',fixsize);
    pO4 = plot(0,0,'ok','MarkerSize',fixsize);
    
    %plot lines between joints
    p2 = plot([0 xAi(i)],[yO2 yAi(i)],'k-');
    p4 = plot([0 xBi(i)],[0 yBi(i)], 'k-');
    p5 = plot([xBi(i) xCi(i)],[yBi(i) yC],'k-');
    
    xlim(xlimval)
    ylim(ylimval)
    xlabel('x-coordinates [m]')
    ylabel('y-coordinate [m]')
    title('system movement animation')
    drawnow
    movievector(i) = getframe;
end

%% write movie
clear myWriter
delete hw7.avi
myWriter = VideoWriter('hw7');
myWriter.FrameRate = 20;
open(myWriter)
writeVideo(myWriter,movievector)
close(myWriter);