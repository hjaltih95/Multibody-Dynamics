clear all; close all; clc;
syms L g m I x1 y1 x2 y2 xd1 yd1 xd2 yd2 xdd phi1 phid1 phi2 phid2 phidd1 phidd2 t lambda rhosym;

% Initializing mass matrix
M = diag([m m I m m I]);

% These are the generalized coordinates q and its derivatives
q = [phi1; phi2];
qd = [phid1; phid2];
qdd = [phidd1; phidd2];


% The coordinates of x1, y1, x2, y2 in generalized coordinates
x1 = L/2*cos(phi1);
y1 = L/2*sin(phi1);
x2 = L*cos(phi1) + L/2*cos(phi2);
y2 = L*sin(phi1) + L/2*sin(phi2);


% qkd is the derivative of qk
qk = [x1; y1; phi1; x2; y2; phi2];
qkd = jacobian(qk,q)*qd;
qkd = simplify(qkd);
qkdd = simplify(jacobian(qkd,q))*qd + simplify(jacobian(qkd,qd))*qdd;


% Constraint equations for the spring
xD = -L/2;
yD = 0;
xE = x1 + (1/6)*cos(phi1)*L;
yE = y1 + (1/6)*sin(phi1)*L;


% Spring stiffness
k = (15/2)*((m*g)/L);
Ls = sqrt((xE - xD)^2 + (yE - yD)^2);
L0 = (2/3)*L;
omega = -(60*2*pi)/60;
dMotor = phi1 - omega*t;
dImpact = x2 + (L/2)*cos(phi2);
Cs = jacobian(dMotor,q);
Cd = Cs*qd;
C2 = simplify(jacobian(Cd, q)*qd);
CsImpact = jacobian(dImpact,q);
CdImpact = CsImpact*qd;


% Potential engergy
V = m*g*(L/2-x1) + m*g*((3*L)/2-x2);
Vspring = V + (1/2)*k*(Ls-L0)^2;


% Kinetic energy
T = sum(1/2*M*qkd.^2);

% Lagrange formula, there are no external forces
Qj = [0; 0];

% Potential energy
Vd = simplify(jacobian(V,q)).';
Vdspring = simplify(jacobian(Vspring,q)).';

% kinetic energy
Td = simplify(jacobian(T,q)).';
Tdd = simplify(jacobian(T,qd)).';
Tddd = simplify(jacobian(Tdd,q)*qd) + simplify(jacobian(Tdd,qd)*qdd);

% EoM
eom = simplify(Tddd - Td + Vd - Qj);

eomSpring = simplify(Tddd - Td + Vdspring - Qj);

eomMotor = simplify(Tddd - Td + Vd - Qj + (Cs*lambda).');

eomImpact = simplify(Tddd - Td + Vd - Qj + (CsImpact*rhosym).');



[LHS, RHS] = equationsToMatrix(eom, qdd);
[LHSSpring, RHSSpring] = equationsToMatrix(eomSpring, qdd);
[LHSMotor, RHSMotor] = equationsToMatrix(eomMotor, [qdd.' lambda]);
[LHSImpact, RHSImpact] = equationsToMatrix(eomImpact, [qdd.' rhosym]);


% Matrix for motor question
MatrixMotor = [LHSMotor; Cs 0];
RHSMotor = [RHSMotor; -C2];


% For different e's
MatrixImpact = [LHSImpact(1:2,1:2) CsImpact.'; CsImpact 0];
RHSImpact0 = [LHSImpact(1:2,1:2)*qd; -0*CsImpact*qd];
RHSImpact05 = [LHSImpact(1:2,1:2)*qd; -0.5*CsImpact*qd];
RHSImpact1 = [LHSImpact(1:2,1:2)*qd; -1*CsImpact*qd];
toSub = simplify(LHS\RHS);
toSubSpring = simplify(LHSSpring\RHSSpring);
toSubMotor = simplify(MatrixMotor\RHSMotor);
toSubImpact0 = simplify(MatrixImpact\RHSImpact0);
toSubImpact05 = simplify(MatrixImpact\RHSImpact05);
toSubImpact1 = simplify(MatrixImpact\RHSImpact1);


% Setting constants
rho = 670;
w = 0.05;
t = 0.003;
L = 0.6;
m = w*t*L*rho;
I = (1/12)*m*(L^2+w^2);
g = 9.81;

Constants = {L, m, I, g};


AnswerC_B = getAns(M, toSub, qkdd, [pi/2 pi/2 0 0], 'a');

AnswerC_C = getAns(M, toSub, qkdd, [0 0 0 0], 'a');

InitialAngVelocityD = (30*2*pi)/60;
AnswerC_D = getAns(M, toSub, qkdd, [0 0 InitialAngVelocityD InitialAngVelocityD], 'a');

AnswerC_D_vel = getAns(M, toSub, qkd, [0 0 InitialAngVelocityD InitialAngVelocityD], 'v');

AnswerD = getAns(M, toSubSpring, qkdd, [pi/2 pi/2 0 0], 'a');

AnswerE = getAns(M, toSubMotor, qkdd, [pi/2 pi/2 omega omega], 'l');
AnswerE = [AnswerE; AnswerE(7)*omega];

AnswerF0 = getAns(M, toSubImpact0, qkd, [pi/2 pi/2 -omega -omega], 'k');
AnswerF05 = getAns(M, toSubImpact05, qkd, [pi/2 pi/2 -omega -omega], 'k');
AnswerF1 = getAns(M, toSubImpact1, qkd, [pi/2 pi/2 -omega -omega], 'k');

initialConditions = [pi/2 pi/2 -omega -omega];
phi1 = initialConditions(1);
phi2 = initialConditions(2);
phid1 = initialConditions(3);
phid2 = initialConditions(4);
TBefore = vpa(subs(sum((1/2)*M*qkd.^2)), 3);
qdd = subs(toSubImpact0);

phid1 = qdd(1);
phid2 = qdd(2);
TAfter0 = vpa(subs(sum((1/2)*M*AnswerF0(1:6,1).^2)), 3);
Tdif0 = double(vpa(TAfter0 - TBefore, 3));

qdd = subs(toSubImpact05);
phid1 = qdd(1);
phid2 = qdd(2);
TAfter05 = vpa(subs(sum((1/2)*M*AnswerF05(1:6,1).^2)), 3);
Tdif05 = double(vpa(TAfter05 - TBefore, 3));

qdd = subs(toSubImpact1);
phid1 = qdd(1);
phid2 = qdd(2);
TAfter1 = vpa(subs(sum((1/2)*M*AnswerF1(1:6,1).^2)), 3);
Tdif1 = double(vpa(TAfter1 - TBefore, 3));


% Displaying a table for C_B, C_C and C_D accelerations
Variables = {'xdd1'; 'ydd1'; 'phidd1'; 'xdd2'; 'ydd2'; 'phidd2'};
Units = {'m/s^2'; 'm/s^2'; 'rad/s^2'; 'm/s^2'; 'm/s^2'; 'rad/s^2'};
Table = table(Variables, AnswerC_B, AnswerC_C, AnswerC_D, Units)
% Displaying a table for C_D velocities
Variables = {'xd1'; 'yd1'; 'phid1'; 'xd2'; 'yd2'; 'phid2'};
Units = {'m/s'; 'm/s'; 'rad/s'; 'm/s'; 'm/s'; 'rad/s'};
Table = table(Variables, AnswerC_D_vel, Units)
% Displaying a table for D accelerations
Variables = {'xdd1'; 'ydd1'; 'phidd1'; 'xdd2'; 'ydd2'; 'phidd2'};
Units = {'m/s^2'; 'm/s^2'; 'rad/s^2'; 'm/s^2'; 'm/s^2'; 'rad/s^2'};
Table = table(Variables, AnswerD, Units)
% Displaying a table for E accelerations
Variables = {'xdd1'; 'ydd1'; 'phidd1'; 'xdd2'; 'ydd2'; 'phidd2'; 'lambda';
'Power'};
Units = {'m/s^2'; 'm/s^2'; 'rad/s^2'; 'm/s^2'; 'm/s^2'; 'rad/s^2'; 'Nm';
'W'};
Table = table(Variables, AnswerE, Units)
% Displaying a table for F velocities
Variables = {'xd1'; 'yd1'; 'phid1'; 'xd2'; 'yd2'; 'phid2'; 'rho'; 'dT';
'W'};
Units = {'m/s'; 'm/s'; 'rad/s'; 'm/s^2'; 'm/s'; 'rad/s'; 'kgm/s'; 'J';
'J'};
Table = table(Variables, [AnswerF0; Tdif0; -Tdif0], [AnswerF05; Tdif05; 
Tdif05], [AnswerF1; Tdif1; -Tdif1], Units);
Table.Properties.VariableNames = {'Variables' 'AnswerF_e_0' 'AnswerF_e_half' 'AnswerF_e_1' 'Units'}
% Super nice function to create the final answers.
function Answer = getAns(M, toSub, conditions, initialConditions, mode)

 phi1 = initialConditions(1);
 phi2 = initialConditions(2);
 phid1 = initialConditions(3);
 phid2 = initialConditions(4);

 if mode == 'l'
 qdd = subs(toSub);
 phidd1 = qdd(1);
 phidd2 = qdd(2);

 Answer = [double(vpa(subs(conditions), 6)); double(vpa(qdd(3),6))];
 end
 if mode == 'k'
 qdd = subs(toSub);
 phid1 = qdd(1);
 phid2 = qdd(2);

 Answer = [double(vpa(subs(conditions), 6)); double(vpa(qdd(3),6))];

 end
 if mode == 'a'
 qdd = subs(toSub);
 phidd1 = qdd(1);
 phidd2 = qdd(2);

 Answer = double(vpa(subs(conditions), 6));
 end
 if mode == 'v'
 Answer = double(vpa(subs(conditions), 6));

 end
end