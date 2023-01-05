%% passive spring
clc;clear;
w = 0.05; h = 0.003; p = 670; g = 9.81; l1 = 0.6; l2 = 0.6;
m = p*w*h*l2;
% I = p*w*h*l2^3/12;
I = m*(l1^2+w^2)/12;
k = (15/2)*(m*g/l1);
l0 = 2*l1/3;

% initial condition
phi1 = pi/2; phi2 = pi/2; phid1 = 0; phid2 = 0; x1 = 0; y1 = l1/2;  

Cs = sqrt((x1+l1*sin(phi1)/6+l1/2)^2 + (y1+l1*cos(phi1)/6)^2) - l0;
sita = k*Cs;
%   define matrix
M = [m 0 0 0 0 0;
     0 m 0 0 0 0;
     0 0 I 0 0 0;
     0 0 0 m 0 0;
     0 0 0 0 m 0;
     0 0 0 0 0 I;
     ];

Fg = [m*g 0 0 m*g 0 0]';

Cx =[1,  0,  (l1*sin(phi1))/2, 0, 0,             0;
 0,  1, -(l1*cos(phi1))/2, 0, 0,                 0;
-1,  0,  (l1*sin(phi1))/2, 1, 0,  (l2*sin(phi2))/2;
 0, -1, -(l1*cos(phi1))/2, 0, 1, -(l2*cos(phi2))/2];

Cl =[ (l1 + 2*x1 - (l1*sin(phi1))/3)/(2*((l1/2 + x1 - (l1*sin(phi1))/6)^2 + (y1 + (l1*cos(phi1))/6)^2)^(1/2));
    (2*y1 + (l1*cos(phi1))/3)/(2*((l1/2 + x1 - (l1*sin(phi1))/6)^2 + (y1 + (l1*cos(phi1))/6)^2)^(1/2));
     -(l1*(l1*cos(phi1) + 2*x1*cos(phi1) + 2*y1*sin(phi1)))/(2*((3*l1 + 6*x1 - l1*sin(phi1))^2 + (6*y1 + l1*cos(phi1))^2)^(1/2));
    0; 0; 0];

C2 = [                     (l1*phid1^2*cos(phi1))/2;
                            (l1*phid1^2*sin(phi1))/2;
 (l1*cos(phi1)*phid1^2)/2 + (l2*cos(phi2)*phid2^2)/2;
 (l1*sin(phi1)*phid1^2)/2 + (l2*sin(phi2)*phid2^2)/2;];

V = [M Cx'; Cx zeros(4)];
b = [Fg - Cl*sita; -C2];

% calculation
X = V\b

%% positive element
w = 0.05; h = 0.003; p = 670; g = 9.81; l1 = 0.6; l2 = 0.6;
m = p*w*h*l2;
% I = p*w*h*l2^3/12;
I = m*(l1^2+w^2)/12;
% Counterclockwise is defined as positive direction
omega = -2*pi;

% initial condition
phi1 = pi/2; phi2 = pi/2; phid1 = omega; phid2 = omega; x1 = 0; y1 = l1/2; x2 = 0; y2 = l1+l2/2;

M = [m 0 0 0 0 0;
     0 m 0 0 0 0;
     0 0 I 0 0 0;
     0 0 0 m 0 0;
     0 0 0 0 m 0;
     0 0 0 0 0 I;
     ];

Fg = [m*g 0 0 m*g 0 0]';

Cx =[1,  0,  (l1*sin(phi1))/2, 0, 0,             0;
 0,  1, -(l1*cos(phi1))/2, 0, 0,                 0;
-1,  0,  (l1*sin(phi1))/2, 1, 0,  (l2*sin(phi2))/2;
 0, -1, -(l1*cos(phi1))/2, 0, 1, -(l2*cos(phi2))/2;
 0, 0, 1, 0, 0, 0];

C2 = [                     (l1*phid1^2*cos(phi1))/2;
                            (l1*phid1^2*sin(phi1))/2;
 (l1*cos(phi1)*phid1^2)/2 + (l2*cos(phi2)*phid2^2)/2;
 (l1*sin(phi1)*phid1^2)/2 + (l2*sin(phi2)*phid2^2)/2;
 0];

V = [M Cx'; Cx zeros(5)];
b = [Fg; -C2];

% calculation
X = V\b


% Power: Torque*angular velocity
P = omega*X(10);
%% impact
clc;clear;
w = 0.05; h = 0.003; p = 670; g = 9.81; l1 = 0.6; l2 = 0.6;
m = p*w*h*l2;
% I = p*w*h*l2^3/12;
I = m*(l1^2+w^2)/12;
omega = 2*pi;

syms phi1 phi2 phi1d phi2d
x1 = cos(phi1)*l1/2;
y1 = sin(phi1)*l1/2;
x2 = cos(phi1)*l1+cos(phi2)*l2/2;
y2 = sin(phi1)*l1+cos(phi2)*l2/2;
xx = [x1;y1;x2;y2];
q = [phi1; phi2];
qd = [phi1d;phi2d];

% velocities calculation, velocities before impact!
xdx = jacobian(xx,q)*qd


% initial condition
phi1 = pi/2; phi2 = pi/2; phid1 = omega; phid2 = omega; x1 = 0; y1 = l1/2; x2 = 0; y2 = l1+l2/2; 




xd1 = -sin(phi1)*l1*phid1/2; yd1 = cos(phi1)*l1*phid1/2; xd2 = -2*l2*phid1*sin(phi1)/2-sin(phi2)*l2*phid2/2; yd2 = 2*l2*phid1*cos(phi1)/2 + cos(phi2)*l2*phid2/2;
x_m = [xd1; yd1; phid1; xd2; yd2; phid2];
e1 = 0;
e2 = 0.5;
e3 = 1;

M = [m 0 0 0 0 0;
     0 m 0 0 0 0;
     0 0 I 0 0 0;
     0 0 0 m 0 0;
     0 0 0 0 m 0;
     0 0 0 0 0 I;
     ];
 
 Fg = [m*g 0 0 m*g 0 0]';
 
Cx =[1,  0,  (l1*sin(phi1))/2, 0, 0,             0;
 0,  1, -(l1*cos(phi1))/2, 0, 0,                 0;
-1,  0,  (l1*sin(phi1))/2, 1, 0,  (l2*sin(phi2))/2;
 0, -1, -(l1*cos(phi1))/2, 0, 1, -(l2*cos(phi2))/2;
 0,  0,                 0, 1, 0, -(l2*sin(phi2))/2];



C2 = [                     (l1*phid1^2*cos(phi1))/2;
                            (l1*phid1^2*sin(phi1))/2;
 (l1*cos(phi1)*phid1^2)/2 + (l2*cos(phi2)*phid2^2)/2;
 (l1*sin(phi1)*phid1^2)/2 + (l2*sin(phi2)*phid2^2)/2;
 -(l2*phid2^2*cos(phi2))/2;
 -(l2*phid2^2*sin(phi2))/2];

V = [M Cx'; Cx zeros(5)];

% three configuration, see 12.88 book!
b1 = [M*x_m; -e1*Cx*x_m];
b2 = [M*x_m; -e2*Cx*x_m];
b3 = [M*x_m; -e3*Cx*x_m];

% calculation
X_1 = V\b1
X_2 = V\b2
X_3 = V\b3