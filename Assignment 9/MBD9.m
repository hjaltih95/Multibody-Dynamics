%% a) Find Rotation Matrix from the body to the inertial frame nRb
% Body frame
e_bx = [0.768;0.024;0.640];
e_by = [-0.424;0.768;0.480];
e_bz = [-0.480; -0.640; 0.600];
e_b = [e_bx, e_by, e_bz];

% Inertial frame
e_x = [1;0;0];
e_y = [0;1;0];
e_z = [0;0;1];
e = [e_x, e_y, e_z];

% Rotation matrix from body to inertial frame
R = e*e_b;

%% b) Determine the Euler angles  (z-x-z)
syms phi theta psi
ang = [phi; theta; psi];

% Find first the rotation marixes about z-x-z axes, i.e. phi,theta and psi,
% individually

R_phi = [cos(phi) -sin(phi) 0;          % Same rotation matrix as in 19.26 in the book
        sin(phi) cos(phi) 0;
        0 0 1];

R_theta = [1 0 0;                       % Same rotation matrix as in 19.28 in the book
           0 cos(theta) -sin(theta);
           0 sin(theta) cos(theta)];
       
R_psi = [cos(psi) -sin(psi) 0;          % Same rotation matrix as in 19.30 in the book
        sin(psi) cos(psi) 0;
        0 0 1];

    
% The final rotation matrix

R_Eul = R_phi*R_theta*R_psi;

matlabFunction(R_Eul,'File','Rhw9');

% Now we find the value of the angles

theta = acos(R(3,3));
phi = asin(R(1,3)/sin(theta));
psi = asin(R(3,1)/sin(theta));

Eul = [phi;theta;psi];

%% c) Finding the angular velocities in the body frame

% Angular velocities in the inertia frame at t=0
w0_i = [7.67952;0.23936;6.40060];

% Calculate angular velocities in body frame, using equation 19.60
w0_b = inv(R)*w0_i;


%% d) Find the rate of change of Euler angles
% we use equation 19.70 to determine the rate of change
B = [sin(psi)*sin(theta) cos(psi) 0;
    cos(psi)*sin(theta) -sin(psi) 0;
    cos(theta) 0 1];

dEul = inv(B)*w0_b;

%% e) Find the inertia matrix bIc in the body fixed frame.

m = 60;    % [kg]

l_x = 0.4;  % [m]
l_y = 1.2;  % [m] 
l_z = 0.3;  % [m]
    
% Moment of inertia around the center of mass is 1/12*ml^2, and we use the
% inertia matrix form in equation 18.28, i.e. we align the principal axes of
% the body such that we get a diagonal mass moment of inertia, i.e. such
% that the object is symmetric around its pricipal axes.

bIx = 1/12*m*(l_y^2+l_z^2);
bIy = 1/12*m*(l_x^2+l_z^2);
bIz = 1/12*m*(l_x^2+l_y^2);

bIc = diag([bIx bIy bIz]);
    

%% f) Equation of motion
syms bwx bwy bwz Ixx Iyy Izz

wb = [bwx; bwy; bwz];
I = diag([Ixx Iyy Izz]);

dwb = inv(I)*(-cross(wb,(I*wb)));

%% g) Numerical Integration with the Runge-Kutta 4th order method

% setup
time=60;
nn=16;
N=2.^nn;
h=time./N;

t= 0;
y0 = [phi; theta; psi; w0_b];
ang = [phi; theta; psi; w0_b];
for j=1:N
    k1 = qa(y0);
    k2=qa(y0 + h/2*k1);
    k3=qa(y0 + h/2*k2);
    k4=qa(y0 + h*k3);
    qn= y0+1/6*h*(k1+2*k2+2*k3+k4);
    
    y0 = qn;
    q_n(j,:)=qn;
    t = t+h;
    T(j)=t;
end
    
%% Plots
q_plot = [ang';q_n];
T_plot = [0;T'];

figure(1)
plot(T_plot,q_plot(:,1),T_plot,q_plot(:,2),T_plot,q_plot(:,3));
xlabel('Time [s]')
ylabel('Euler angles [rad]')
legend('Phi [\phi]', 'Theta [\theta]', 'Psi [\psi]')
title('Euler angles vs time')

figure(2)
plot(T_plot,q_plot(:,4),T_plot,q_plot(:,5),T_plot,q_plot(:,6));
xlabel('Time [s]')
ylabel('Angular velocities in the body frame [rad/s]')
legend('Omgea_x [\omega_x]', 'Omega_y [\omega_y]', 'Omega_z [\omega_z]')
title('Angular velocities in the body frame')

figure(3)
subplot(311)
plot(T_plot,q_plot(:,1))
title('Euler angles vs time')
ylabel('Euler angles [rad]')
legend('Phi [\phi]')

subplot(312)
plot(T_plot,q_plot(:,2))
ylabel('Euler angles [rad]')
legend('Theta [\theta]')

subplot(313)
plot(T_plot,q_plot(:,3))
ylabel('Euler angles [rad]')
legend('Psi [\psi]')
xlabel('Time [s]')

figure(4)
subplot(311)
plot(T_plot,q_plot(:,4))
title('Angular velocities in the body frame')
ylabel('Omega_x (\omega_x) in the body frame [rad/s]')
legend('Omgea_x [\omega_x]')

subplot(312)
plot(T_plot,q_plot(:,5))
ylabel('Omega_y (\omega_y) in the body frame [rad/s]')
legend( 'Omega_y [\omega_y]')

subplot(313)
plot(T_plot,q_plot(:,6))
ylabel('Omega_z (\omega_z) in the body frame [rad/s]')
legend('Omega_z [\omega_z]')
xlabel('Time [s]')

figure(9)
plot(T_plot,q_plot(:,2),T_plot,q_plot(:,3))
xlabel('Time [s]')
ylabel('Angular velocities in the body frame [rad/s]')
legend('Omega_y [\omega_y]', 'Omega_z [\omega_z]')
title('Angular velocities in the body frame')
%% h) 3D trajectory of point p
p = [(l_x)/2;0;0];

p_I = R_Euler(p,q_plot(:,1:3));

figure(5)
plot3(p_I(:,1),p_I(:,2),p_I(:,3))
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Trajectory of point p in the inertial frame')

%% i) plot the angular velocities in the inertial frame

w_i = R_Euler_vel(q_plot);

figure(6)
plot(T_plot,w_i(:,1),T_plot,w_i(:,2),T_plot,w_i(:,3))
title('Angular velocities in the inertial frame')
ylabel('Angular velocities[rad/s]')
xlabel('Time [s]')
legend('Omgea_x [\omega_x]', 'Omega_y [\omega_y]', 'Omega_z [\omega_z]')


figure(7)
subplot(311)
plot(T_plot,w_i(:,1))
title('Angular velocities in the inertial frame')
ylabel('Omega_x, \omega_x in the inertial frame [rad/s]')
legend('Omgea_x [\omega_x]')

subplot(312)
plot(T_plot,w_i(:,2))
ylabel('Omega_y, \omega_y in the inertial frame  [rad/s]')
legend('Omgea_y [\omega_y]')

subplot(313)
plot(T_plot,w_i(:,3))
ylabel('Omega_z, \omega_z in the inertial frame  [rad/s]')
legend('Omgea_z [\omega_z]')
xlabel('Time [s]')

%% i) Angular momentum vector in the inertial frame

momentum_i = R_Euler_mom(q_plot);

figure(8)
plot(T_plot,momentum_i(:,1),T_plot,momentum_i(:,2),T_plot,momentum_i(:,3))
ylabel('Angular momentum [kgm^2/s]')
xlabel('Time [s]')
title('Angular momentum expressed in the inertial frame as a function of time')
legend('H_x', 'H_y', 'H_z')

%% j) Invariants

% distance between two points on the satellite.
p1 = [0;l_y;0];
p2 = [l_x;0;0];

P_1 = R_Euler(p1,q_plot(:,1:3));
P_2 = R_Euler(p2,q_plot(:,1:3));

for i = 1:length(T_plot)
    dist(i,:) = norm(P_1(i,:)-P_2(i,:));
end

figure(10)
plot(T_plot,dist)
ylabel('Distance [m]')
xlabel('Time [s]')
title('Distance between 2 points on the satellite')
yaxis([0 1.5])

%% Kinetic energy of the system

T = R_Euler_KE(q_plot);

figure(11)
plot(T_plot,T)
xlabel('Time [s]')
ylabel('Energy of the system')
title('Kineric energy of the system for each timestep')
%yaxis([0 850])

%% Functions

function state_d = qa(y)
phi = y(1);
theta = y(2);
psi = y(3);

wb = y(4:6);

% Find the rate of change in angular velocities
m = 60;    % [kg]

l_x = 0.4;  % [m]
l_y = 1.2;  % [m] 
l_z = 0.3;  % [m]


bIx = 1/12*m*(l_y^2+l_z^2);
bIy = 1/12*m*(l_x^2+l_z^2);
bIz = 1/12*m*(l_x^2+l_y^2);

bIc = diag([bIx bIy bIz]);

dwb = inv(bIc)*(-cross(wb,(bIc*wb)));

% Find the rate of change for Euler angles
B = [sin(psi)*sin(theta) cos(psi) 0;
    cos(psi)*sin(theta) -sin(psi) 0;
    cos(theta) 0 1];

dEul = inv(B)*wb;

state_d = [dEul;dwb];
end



function position = R_Euler(p,angles)

for i = 1:length(angles)
    phi = angles(i,1);
    theta = angles(i,2);
    psi = angles(i,3);
    
    R_Eul = Rhw9(phi,psi,theta);
    
    position(i,:) = R_Eul*p;
end

end

function angvel = R_Euler_vel(angles)

for i = 1:length(angles)
    phi = angles(i,1);
    theta = angles(i,2);
    psi = angles(i,3);
    
    wx = angles(i,4);
    wy = angles(i,5);
    wz = angles(i,6);
    
    wb = [wx;wy;wz];
    
    R_Eul = Rhw9(phi,psi,theta);
    
    angvel(i,:) = R_Eul*wb;
end

end

function momentum = R_Euler_mom(angles)

m = 60;    % [kg]

l_x = 0.4;  % [m]
l_y = 1.2;  % [m] 
l_z = 0.3;  % [m]


bIx = 1/12*m*(l_y^2+l_z^2);
bIy = 1/12*m*(l_x^2+l_z^2);
bIz = 1/12*m*(l_x^2+l_y^2);

bIc = diag([bIx bIy bIz]);


for i = 1:length(angles)
    phi = angles(i,1);
    theta = angles(i,2);
    psi = angles(i,3);
    
    wx = angles(i,4);
    wy = angles(i,5);
    wz = angles(i,6);
    
    wb = [wx;wy;wz];
    
    R_Eul = Rhw9(phi,psi,theta);
    
    momentum(i,:) = R_Eul*bIc*wb;
end

end

function kenergy = R_Euler_KE(angles)

m = 60;    % [kg]

l_x = 0.4;  % [m]
l_y = 1.2;  % [m] 
l_z = 0.3;  % [m]


bIx = 1/12*m*(l_y^2+l_z^2);
bIy = 1/12*m*(l_x^2+l_z^2);
bIz = 1/12*m*(l_x^2+l_y^2);

bIc = diag([bIx bIy bIz]);

for i = 1:length(angles)
    
    wx = angles(i,4);
    wy = angles(i,5);
    wz = angles(i,6);
    
    wb = [wx;wy;wz];
    
    kenergy(i,:) = 1/2*wb'*bIc*wb;
end

end
