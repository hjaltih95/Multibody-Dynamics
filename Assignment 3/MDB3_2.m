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
ks=(15/2)*(m1*g/l1)
l0=2/3*l1
%% Initial position and velocities
phi1=pi/2;  
phi2=pi/2;  
phid1=-60/60*2*pi;    
phid2=-60/60*2*pi;
x1=0;
y1=1/2*l1;
omega=-60/60*2*pi
%% Set up DAEs a)
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
    
ls=sqrt((x1+1/6*l1*cos(phi1)+l1/2)^2+(y1+1/6*l1*sin(phi1)-0)^2);

C_si = (1/ls)*[x1+l1/6*cos(phi1)+l1/2;
        y1+l1/6*sin(phi1);
        (-1/12*l1)*(l1*sin(phi1)-2*y1*cos(phi1)+2*x1*sin(phi1));
        0; 0; 0];

C_s=sqrt((x1+1/6*l1*cos(phi1)+l1/2)^2+(y1+1/6*l1*sin(phi1)-0)^2)-l0;

F_i=[m1*g 0 0 m2*g 0 0]';
Mat=[M C_ki; C_kj zeros(4)];
Y=[F_i-C_si*ks*C_s; -C_kij];
Sol=inv(Mat)*Y
xddot=[Sol(1:6)]
lambda=[Sol(7:10)]

%% DAEs for question b)

m=[m1 m1 I1 m2 m2 I2];
M=diag(m);

C_kj=[1 0 (l1*sin(phi1))/2 0 0 0;
      0 1 -(l1*cos(phi1))/2 0 0 0;
      -1 0 (l1*sin(phi1))/2 1 0 (l2*sin(phi2))/2;
      0 -1 -(l1*cos(phi1))/2 0 1 -(l2*cos(phi2))/2;
      0 0 1 0 0 0];
  
C_ki=C_kj';

C_kij=[(l1*phid1^2*cos(phi1))/2;
       (l1*phid1^2*sin(phi1))/2;
        (l1*cos(phi1)*phid1^2)/2 + (l2*cos(phi2)*phid2^2)/2;
        (l1*sin(phi1)*phid1^2)/2 + (l2*sin(phi2)*phid2^2)/2;
        0];

F_i=[m1*g 0 0 m2*g 0 0]';
Mat=[M C_ki; C_kj zeros(5)];
Y=[F_i; -C_kij];
Sol=inv(Mat)*Y
xddot=[Sol(1:6)]
lambda=[Sol(7:11)]

%% DAEs for c)

m=[m1 m1 I1 m2 m2 I2];
M=diag(m);

C_kj=[1 0 (l1*sin(phi1))/2 0 0 0;
      0 1 -(l1*cos(phi1))/2 0 0 0;
      -1 0 (l1*sin(phi1))/2 1 0 (l2*sin(phi2))/2;
      0 -1 -(l1*cos(phi1))/2 0 1 -(l2*cos(phi2))/2;
      0 0 1 0 0 0];
  
C_ki=C_kj';

C_kij=[(l1*phid1^2*cos(phi1))/2;
       (l1*phid1^2*sin(phi1))/2;
        (l1*cos(phi1)*phid1^2)/2 + (l2*cos(phi2)*phid2^2)/2;
        (l1*sin(phi1)*phid1^2)/2 + (l2*sin(phi2)*phid2^2)/2;
        0];

F_i=[m1*g 0 0 m2*g 0 0]';
Mat=[M C_ki; C_kj zeros(5)];
Y=[F_i; -C_kij];
Sol=inv(Mat)*Y
xddot=[Sol(1:6)]
lambda=[Sol(7:11)]