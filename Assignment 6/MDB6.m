%% Set up EoM
syms alpha beta
syms alphadot betadot 

q=[alpha; beta];
qd=[alphadot; betadot];

rho=670;
l=0.6;
V=l*0.05*0.003;
m=rho*V;
I=1/12*m*l^2;
g=9.81;
Qa=0;
Qb=0;

M=[I+5/4*m*l^2, 1/2*m*l^2*cos(alpha-beta);
    1/2*m*l^2*cos(alpha-beta), I+1/4*m*l^2];
F=[Qa-3/2*m*g*l*sin(alpha)-1/2*m*l^2*betadot^2*sin(alpha-beta);
    Qb-1/2*m*g*l*sin(beta)+1/2*m*l^2*alphadot^2*sin(alpha-beta)];


qdd=M\F;

matlabFunction(qdd,'File','qddhw6')

%% Setup
T=3;
nn=5:1:20;
N=2.^nn;
h=T./N;

format('Long')
%% Euler method

% Integration of all time steps

for i=1:length(nn)
    
    y0=[pi/2;pi/2;0;0];
    alpha=y0(1);
    beta=y0(2);
    alphadot=y0(3);
    betadot=y0(4);
    
    for j = 1:N(i)
    
        p=[alpha; beta];
        v=[alphadot; betadot];
    
        pn=p+h(i)*v;
        vn=v+h(i)*qddhw6(alpha,alphadot,betadot,beta);
    
        alpha = pn(1);
        beta = pn(2);
        alphadot = vn(1);
        betadot = vn(2);
    
    end
    q_end_h(i,:) = [alpha, beta];
end

% Global error 
ph=1;
C2=10^(-16);

for i=1:length(nn)-1
    Da_eul(i) = q_end_h(i+1,1)-q_end_h(i,1);
    Db_eul(i) = q_end_h(i+1,2)-q_end_h(i,2);
    
    Eta_euler(i)=abs(1/(2^ph-1)*Da_eul(i));
    Etb_euler(i)=abs(1/(2^ph-1)*Db_eul(i));
    
    Eround(i)=abs(C2/h(i));
    
    E_tota_euler(i)=Eta_euler(i)+Eround(i);
    E_totb_euler(i)=Etb_euler(i)+Eround(i);
    
end

%% Heun method

% Integration of all time steps
for i=1:length(nn)
    
    alpha=pi/2;
    beta=pi/2;
    alphadot=0;
    betadot=0;
    
    for j = 1:N(i)
        
        y0=[alpha;beta;alphadot;betadot];
        
        ystar= y0+h(i)*qa(y0);
        qn=y0+h(i)/2*(qa(y0)+qa(ystar));
        
        alpha = qn(1);
        beta = qn(2);
        alphadot = qn(3);
        betadot = qn(4);
    
    end
    q_end_h(i,:) = [alpha, beta];
end

% Global error 
ph=2;
C2=10^(-16);

for i=1:length(nn)-1
    Da_heun(i) = q_end_h(i+1,1)-q_end_h(i,1);
    Db_heun(i) = q_end_h(i+1,2)-q_end_h(i,2);
    
    Eta_heun(i)=abs(1/(2^ph-1)*Da_heun(i));
    Etb_heun(i)=abs(1/(2^ph-1)*Db_heun(i));
    
    Eround(i)=abs(C2/h(i));
    
    E_tota_heun(i)=Eta_heun(i)+Eround(i);
    E_totb_heun(i)=Etb_heun(i)+Eround(i);
    
end


%% Runge-Kutta 3rd order

% Intergration of all time steps
for i=1:length(nn)
    
    alpha=pi/2;
    beta=pi/2;
    alphadot=0;
    betadot=0;
    
    for j=1:N(i)
        y0=[alpha;beta;alphadot;betadot];
        k1=qa(y0);
        k2=qa(y0 + h(i)/2*k1);
        k3=qa(y0 + h(i)*(2*k2-k1));
        
        qn= y0+1/6*h(i)*(k1+4*k2+k3);
        
        alpha=qn(1);
        beta=qn(2);
        alphadot=qn(3);
        betadot=qn(4);
    end
    q_end_rk3(i,:)=[alpha,beta];
end

% Global error 
p_rk3=3;
C2=10^(-16);

for i=1:length(nn)-1
    Da_rk3(i) = q_end_rk3(i+1,1)-q_end_rk3(i,1);
    Db_rk3(i) = q_end_rk3(i+1,2)-q_end_rk3(i,2);
    
    Eta_rk3(i)=abs(1/(2^p_rk3-1)*Da_rk3(i));
    Etb_rk3(i)=abs(1/(2^p_rk3-1)*Db_rk3(i));
    
    Eround(i)=abs(C2/h(i));
    
    E_tota_rk3(i)=Eta_rk3(i)+Eround(i);
    E_totb_rk3(i)=Etb_rk3(i)+Eround(i);
    
end

%% Runge-Kutta 4th order
% Intergration of all time steps
for i=1:length(nn)
    
    alpha=pi/2;
    beta=pi/2;
    alphadot=0;
    betadot=0;
    
    for j=1:N(i)
        y0=[alpha;beta;alphadot;betadot];
        k1=qa(y0);
        k2=qa(y0 + h(i)/2*k1);
        k3=qa(y0 + h(i)/2*k2);
        k4=qa(y0 + h(i)*k3);
        
        qn= y0+1/6*h(i)*(k1+2*k2+2*k3+k4);
        
        alpha=qn(1);
        beta=qn(2);
        alphadot=qn(3);
        betadot=qn(4);
    end
    q_end_rk4(i,:)=[alpha,beta];
end

% Global error 
p_rk4=4;
C2=10^(-16);

for i=1:length(nn)-1
    Da_rk4(i) = q_end_rk4(i+1,1)-q_end_rk4(i,1);
    Db_rk4(i) = q_end_rk4(i+1,2)-q_end_rk4(i,2);
    
    Eta_rk4(i)=abs(1/(2^p_rk4-1)*Da_rk4(i));
    Etb_rk4(i)=abs(1/(2^p_rk4-1)*Db_rk4(i));
    
    Eround(i)=abs(C2/h(i));
    
    E_tota_rk4(i)=Eta_rk4(i)+Eround(i);
    E_totb_rk4(i)=Etb_rk4(i)+Eround(i);
    
end

%% Plot the error vs step size for the 4 different methods

figure(1)
plot(log10(h(2:end)),log10(E_tota_euler),log10(h(2:end)),log10(E_tota_heun),log10(h(2:end)),log10(E_tota_rk3),log10(h(2:end)),log10(E_tota_rk4), 'linewidth', 1.5)
title('Error with regard to step size for alpha')
xlabel('log_{10}(h)')
ylabel('log_{10}(Error)')
grid on
hold on
plot([log10(10^-6),log10(10^-1)],[log10(10^-6),log10(10^-6)],'k--', 'linewidth', 1);
legend('Euler', 'Heun', 'Runge-Kutta 3rd order', 'Runge-Kutta 4th order', 'Maximum error')



figure(2)
plot(log10(h(2:end)),log10(E_totb_euler), log10(h(2:end)),log10(E_totb_heun), log10(h(2:end)),log10(E_totb_rk3), log10(h(2:end)),log10(E_totb_rk4), 'linewidth', 1.5)
title('Error with regard to step size for beta')
xlabel('log_{10}(h)')
ylabel('log_{10}(Error)')
grid on
hold on
plot([log10(10^-6),log10(10^-1)],[log10(10^-6),log10(10^-6)], 'k--', 'linewidth', 1);
legend('Euler', 'Heun', 'Runge-Kutta 3rd order', 'Runge-Kutta 4th order', 'Maximum error')



%% H max
% maximum allowed error
maxerror=10^(-6);

%Find where the maximum error intercepts with the error function, and
%return the step size at that location
h_alpha_eul=interp1(E_tota_euler(2:end),h(3:end),maxerror);
h_beta_eul=interp1(E_totb_euler(2:end),h(3:end),maxerror);
h_alpha_heun=interp1(E_tota_heun,h(2:end),maxerror);
h_beta_heun=interp1(E_totb_heun,h(2:end),maxerror);
h_alpha_rk3=interp1(E_tota_rk3,h(2:end),maxerror);
h_beta_rk3=interp1(E_totb_rk3,h(2:end),maxerror);
h_alpha_rk4=interp1(E_tota_rk4,h(2:end),maxerror);
h_beta_rk4=interp1(E_totb_rk4,h(2:end),maxerror);

%% Find P for the 4 different methods

% Euler
Y_eul=diff(log10(E_tota_euler));
X_eul=diff(log10(h(2:end)));
avg_slope_euler=nanmean(Y_eul./X_eul);

% Heun
Y_heun=diff(log10(E_tota_heun));
X_heun=diff(log10(h(2:end)));
I_heun=find(Y_heun>0);
Y_heun(I_heun) = [];
X_heun(I_heun) = [];
avg_slope_heun=nanmean(Y_heun/X_heun);

% Runge-Kutta 3rd order
Y_rk3=diff(log10(E_tota_rk3));
X_rk3=diff(log10(h(2:end)));
I_rk3=find(Y_rk3>0);
Y_rk3(I_rk3) = [];
X_rk3(I_rk3) = [];
avg_slope_rk3=nanmean(Y_rk3/X_rk3);

% Runge-Kutta 4th order
Y_rk4=diff(log10(E_tota_rk4));
X_rk4=diff(log10(h(2:end)));
I_rk4=find(Y_rk4>0);
Y_rk4(I_rk4) = [];
X_rk4(I_rk4) = [];
avg_slope_rk4=nanmean(Y_rk4/X_rk4);


%% ODE
tspan = [0 T];
y0 = [pi/2; pi/2; 0; 0];
options = odeset('RelTol', 1e-16, 'AbsTol', 1e-6, 'Stats', 'on');
format long;
odefun = @qa2;
% ode23
[t23,y23] = ode23(odefun,tspan,y0,options);
alpha23=y23(end,1:2);
h23=T/length(t23);

% ode45
[t45,y45] = ode45(odefun,tspan,y0,options);
alpha45=y45(end,1:2);
h45=T/length(t45);

% ode113
[t113,y113] = ode113(odefun,tspan,y0,options);
alpha113=y113(end,1:2);
h113=T/length(t113);

%% Standard First-Order Form
function acc = qa(y)
alpha=y(1);
beta=y(2);
alphadot=y(3);
betadot=y(4);

acc=qddhw6(alpha,alphadot,betadot,beta);

acc = [y(3); y(4); acc];

end
% Function for ODE solver
function acc2 = qa2(t,y)
alpha=y(1);
beta=y(2);
alphadot=y(3);
betadot=y(4);

acc2=qddhw6(alpha,alphadot,betadot,beta);

acc2 = [y(3); y(4); acc2];

end