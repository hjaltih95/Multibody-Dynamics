function [t,y] = rk4step(odefun,t0,y0,h)
t=t0;
y=y0;
k1=feval(odefun,t,y);
k2=feval(odefun,t+h/2,y+(h/2).*k1);
k3=feval(odefun,t+h/2,y+(h/2).*k2);
k4=feval(odefun,t+h,y+h.*k3);
t=t+h;
y=y+(h/6).*(k1+2.*(k2+k3)+k4);
end