%% Initialize variables
rho=670;
l1=0.6;
l2=l1;
V=l1*0.05*0.003;
m1=rho*V;
m2=m1;
I1=1/12*m1*l1^2;
I2=1/12*m2*l2^2;
g=9.81;
%% Initial position and velocities
theta1=0;            % b) pi/2, c) 0, d) 0 
theta2=0;            % b) pi/2, c) 0, d) 0 
thetadot1=pi;        % b) 0, c) 0, d) pi 
thetadot2=pi;        % b) 0, c) 0, d) pi 
%% Set up DAEs
m=[m1 m1 I1 m2 m2 I2];
M=diag(m);

A=[-1 0 1 0; 
   0 -1 0 1; 
   -1/2*l1*sin(theta1) 1/2*l1*cos(theta1) -1/2*l1*sin(theta1) 1/2*l1*cos(theta1);
   0 0 -1 0;
   0 0 0 -1;
   0 0 -1/2*l2*sin(theta2) 1/2*l2*cos(theta2)];

B=[-1 0 -1/2*l1*sin(theta1) 0 0 0;
   0 -1 1/2*l1*cos(theta1) 0 0 0;
   1 0 -1/2*l1*sin(theta1) -1 0 -1/2*l2*sin(theta2);
   0 1 1/2*l1*cos(theta1) 0 -1 1/2*l2*cos(theta2)];


Fg=[m1*g 0 0 m2*g 0 0]'

a=[1/2*l1*thetadot1^2*cos(theta1);
   1/2*l1*thetadot1^2*sin(theta1);
   1/2*l1*thetadot1^2*cos(theta1) + 1/2*l2*thetadot2^2*cos(theta2);
   1/2*l1*thetadot1^2*sin(theta1) + 1/2*l2*thetadot2^2*sin(theta2)];


Mat=[M A; B zeros(4)];
Y=[Fg; a];

Sol=inv(Mat)*Y

xddot=[Sol(1:6)]
Fc=[Sol(7:10)]

%% Part d, calculating velocity
H=[-1/2*l1*sin(theta1)*thetadot1;
    1/2*l1*cos(theta1)*thetadot1;
    1/2*l1*sin(theta1)*thetadot1 + 1/2*l1*sin(theta2)*thetadot2;
    -1/2*l1*cos(theta1)*thetadot1 - 1/2*l1*cos(theta2)*thetadot2];

K=[1 0 0 0;
   0 1 0 0;
   1 0 -1 0;
   0 1 0 -1];

Vel=inv(K)*H