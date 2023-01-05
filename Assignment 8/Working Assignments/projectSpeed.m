function ycorr = projectSpeed(t,ypred)
% project speed on constraint surface, EzyRoller

[conv,sc,Ds,h] = constraint(ypred);
ycorr= ypred;

%Fulfill holonomic constraint

epsdot=Ds*ycorr(7:12);
dx= -Ds.'*((Ds*Ds.')\epsdot);
ycorr(7:12)=ycorr(7:12)+dx;

%Fulfill nonholonomic constraint

dydot = -Ds.'*((Ds*Ds.')\sc);
ycorr(7:12) = vpa(ypred(7:12)+dydot);

end