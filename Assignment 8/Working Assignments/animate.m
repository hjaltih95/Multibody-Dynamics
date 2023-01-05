function animate(T,Y)
% animate motion of chaplygin sleigh
% wb1413, Multibody Dynamics B, Spring 2016.
% Arend L. Schwab 23-Apr-2016
% Copyright (c) 2016 by TU-Delft, the Netherlands.
% Modified by Erdi Akyuz
global A_ M1_ M2_ I1_ I2_ F_ FC_ CD_ B_ C_ D_
X1 = Y(:,1) ; Y1 = Y(:,2) ; P1 = Y(:,3);
X2 = Y(:,4); Y2 = Y(:,5); P2 = Y(:,6);
DX1 = Y(:,7) ; DY1 = Y(:,8) ; DP1 = Y(:,9);
DX2 = Y(:,10); DY2 = Y(:,5); DP2 = Y(:,6);
XA = X1-A_ *cos(P1) ; YA = Y1-A_ *sin(P1);
XB1 = X1+B_ *cos(P1) ; YB1 = Y1+B_ *sin(P1);
XB2 = X2-D_ *cos(P2) ; YB2 = Y2-D_ *sin(P2);
XC = X2+C_ *cos(P2) ; YC = Y2+C_ *sin(P2);
figure
plot(XA,YA)
hold on
plot(XC,YC)
grid on
xlabel('x-position [m]')
ylabel('y-position [m]')
legend('Point A','Point C')
title('Animation EzyRoller')
% axis([min(X1)-A max(X1)+A min(Y1)-A max(Y1)+A ]);
% axis equal
l = plot([XA(1) XB1(1)],[YA(1) YB1(1)]);
k = plot([XB2(1) XC(1)],[YB2(1) YC(1)]);
set(l,'LineWidth',5);
set(l,'Color','K')
set(k,'LineWidth',5);
set(k,'Color','R')
nstep = length(T);
nskip = 1;
for istep = 2:nskip:nstep
set(l,'XData',[XA(istep) XB1(istep)])
set(l,'YData',[YA(istep) YB1(istep)])
set(k,'XData',[XB2(istep) XC(istep)])
set(k,'YData',[YB2(istep) YC(istep)])
drawnow
end
end