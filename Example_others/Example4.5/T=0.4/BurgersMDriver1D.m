% Driver script for solving the 1D Burgers equation using a monotone scheme
clear all
close all
clc
% Set problem parameters
L = 1; FinalTime = 0.4; N = 500; h = L/N;

% Define domain and initial conditions
x = [0:h:1]'; 
u_rk = sin(2*pi*x);
u_rkl = sin(2*pi*x);
u_prk = sin(2*pi*x);
u_prrk = sin(2*pi*x);

%% RK(5,4) scheme with forward Euler time step
k = 0.002;
kappa = 0;
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u_rk] = BurgersM1D(x,u_rk,h,k,kappa);
  time = time+k;
  tstep = tstep+1;
end

%% RK(5,4) scheme with large time step
k = 0.005395;
kappa = 0;
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u_rkl] = BurgersM1D(x,u_rkl,h,k,kappa);
  time = time+k;
  tstep = tstep+1;
end

%% pRK(5,4) scheme
k = 0.01;
kappa = 1/h-1.508/k;
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u_prk] = BurgersM1D(x,u_prk,h,k,kappa);
  time = time+k;
  tstep = tstep+1;
end

%% pRRK(5,4) scheme
k = 0.01;
kappa = 1/h-1.508/k;
c_s = (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+5*64/14293*(k*kappa)^4)/...
      (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+1/24*(k*kappa)^4+64/14293*(k*kappa)^5);
k_hat = c_s*k;
% recalculate timestep of relaxation
[k_res,kappa_res] = relax(h,k_hat,FinalTime);
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  if (time<=fix(FinalTime/k_hat)*k_hat)
%   different "k" corresponding to fix(FinalTime/k_hat) "¡À1"
  [u_prrk] = BurgersM1D(x,u_prrk,h,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u_prrk] = BurgersM1D(x,u_prrk,h,k_res,kappa_res);
  break
  end
end

%% exact solution
for i = 1:N+1
    x(i) = 0 + (i-1)*h;
    xn = bisection(x(i),FinalTime);
    ue(i) = sin(2*pi*xn);
end

%% plots
figure(1)
plot(x,ue,'-',x,u_rk,'.-',x,u_rkl,':',x,u_prk,'--',x,u_prrk,'.:');
hold on;
axis([0,1,-1.2,1.2]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[-1.2:0.4:1.2]);
axis square;
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$u$','interpreter','latex','FontSize',16);
legend('Exact','RK(5,4),\ $\Delta t=0.002$','RK(5,4),\ $\Delta t=0.005395$','pRK(5,4),\ $\Delta t=0.01$','pRRK(5,4),\ $\Delta\hat{t}\approx0.008019$','location','Northeast');
set(gcf,'Units','centimeters','Position',[10 1 15.6 15],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);

figure(2)
plot(x,ue,'-',x,u_rk,'.-',x,u_rkl,':',x,u_prk,'--',x,u_prrk,'.:');
hold on;
axis([0.37,0.53,0.7,1.0]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[-1.2:0.8:1.2]);
axis square;
set (gca,'FontSize',12,'fontweight','demi');
set(gcf,'Units','centimeters','Position',[10 1 15.6 15],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);
