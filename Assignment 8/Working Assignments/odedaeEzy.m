function dyvec = odedaeEzy(t,yvec)

global A_ M1_ M2_ I1_ I2_ F_ FC_ CD_ B_ C_ D_

% setup DAE

M = diag([M1_ M1_ I1_ M2_ M2_ I2_ ]);
[conv,sc,Ds,h] = constraint(yvec);
DAE = [M Ds.';Ds zeros(4,4)];

% determine the applied forces

f = AppliedForces(t,yvec);
rhs = [f; -h; -conv];
lhs = DAE\rhs;

% copy velocities and accelerations into dy

dyvec = yvec;
dyvec(1:6) = yvec(7:12);
dyvec(7:12) = lhs(1:6);

% copy constraint force in global FC

FC_ = [lhs(7);lhs(8);lhs(9);lhs(10)];
end