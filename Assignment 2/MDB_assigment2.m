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
phi1=0;  % b) pi/2, c) 0, d) 0, e) pi/2, f) pi/2, g) 0
phi2=0;  % b) pi/2, c) 0, d) 0, e) pi/2, f)pi/2, g) 0
phid1=0;    % b) 0, c) 0, d) pi, e) 0, f) -pi, g) 0
phid2=0;    % b) 0, c) 0, d) pi, e) 0, f) 0, g) 0 
%% Set up DAEs  b), c) and d)
m=[m1 m1 I1 m2 m2 I2];
M=diag(m);
C_kj=[1 0 (l1*sin(phi1))/2 0 0 0;
      0 1 -(l1*cos(phi1))/2 0 0 0;
      -1 0 (l1*sin(phi1))/2 1 0 (l2*sin(phi2))/2;
      0 -1 -(l1*cos(phi1))/2 0 1 -(l2*cos(phi2))/2];
C_ki=C_kj';
C_kij=[(l1*phid1^2*cos(phi1))/2;
       (l1*phid1^2*sin(phi1))/2;
        (l1*cos(phi1)*phid1^2)/2 + (l2*cos(phi2)*phid2^2)/2;
        (l1*sin(phi1)*phid1^2)/2 + (l2*sin(phi2)*phid2^2)/2];
F_i=[m1*g 0 0 m2*g 0 0]';
Mat=[M C_ki; C_kj zeros(4)];
Y=[F_i; -C_kij];
Sol=inv(Mat)*Y
xddot=[Sol(1:6)]
lambda=[Sol(7:10)]

%% Set up DAEs e), f) and g)
m=[m1 m1 I1 m2 m2 I2];
M=diag(m);

C_kj=[1 0 (l1*sin(phi1))/2 0 0 0;
      0 1 -(l1*cos(phi1))/2 0 0 0;
      -1 0 (l1*sin(phi1))/2 1 0 (l2*sin(phi2))/2;
      0 -1 -(l1*cos(phi1))/2 0 1 -(l2*cos(phi2))/2;
      0 0 0 1 0 -(l2*sin(phi2))/2];
  
C_ki=C_kj';

C_kij=[(l1*phid1^2*cos(phi1))/2;
       (l1*phid1^2*sin(phi1))/2;
        (l1*cos(phi1)*phid1^2)/2 + (l2*cos(phi2)*phid2^2)/2;
        (l1*sin(phi1)*phid1^2)/2 + (l2*sin(phi2)*phid2^2)/2;
        -(l2*phid2^2*cos(phi2))/2];
    
F_i=[m1*g 10 10*l1/2 m2*g 0 0]';    % For g) we put F_i=[m1*g 10 10*l1/2 m2*g 0 0]'
Mat=[M C_ki; C_kj zeros(5)];
Y=[F_i; -C_kij];
Sol=inv(Mat)*Y
xddot=[Sol(1:6)]
lambda=[Sol(7:11)]


