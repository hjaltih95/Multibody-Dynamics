function ycorr = projectSpeed(t,ypred)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
[s,Ds,h] = constraint(ypred);
dydot = -Ds.'*((Ds*Ds.')\s); % [dx1 dy1 dp1 dx2 dy2 dp2]
ycorr = ypred;
ycorr(7:12) = ypred(7:12) + dydot(1:6);
end

