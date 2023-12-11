% Driver script for solving the 1D wave equation using a monotone scheme
clear,clc,close all

%% Set problem parameters
L = 2;              % interval length
FinalTime = 2.0;    % final time
Nvec = 1000;        % spatial grid number
h = L/Nvec;         % spatial step size

%% Define domain and initial conditions
x = [-1:h:1]';
[u_rk] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
[u_rkl] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
[u_prk] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
[u_prrk] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));

%% Compute exact solution
[ue] = waveteste(x-FinalTime,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));

%% RK(5,4) scheme with forward Euler step
k = 0.002;          % time step
kappa = 0;          % stabilization parameter
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u_rk] = LinwaveM1D(x,u_rk,h,k,kappa);
  time = time+k;
  tstep = tstep+1;
end

%% RK(5,4) scheme with large time step
k = 0.004385;
kappa = 0;
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u_rkl] = LinwaveM1D(x,u_rkl,h,k,kappa);
  time = time+k;
  tstep = tstep+1;
end

%% pRK(5,4) scheme
k = 0.5;
kappa = 1/h-1.508/k;
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u_prk] = LinwaveM1D(x,u_prk,h,k,kappa);
  time = time+k;
  tstep = tstep+1;
end

%% pRRK(5,4) scheme
k = 0.5;
kappa = 1/h-1.508/k;
c_s = (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+5*64/14293*(k*kappa)^4)/...
      (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+1/24*(k*kappa)^4+64/14293*(k*kappa)^5);
k_hat = c_s*k;
% recalculate relaxation time step
[k_res,kappa_res] = relax(h,k_hat,FinalTime);
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  if (time<=fix(FinalTime/k_hat)*k_hat)
  [u_prrk] = LinwaveM1D(x,u_prrk,h,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u_prrk] = LinwaveM1D(x,u_prrk,h,k_res,kappa_res);
  break
  end
end

%% figure output
figure(1)
plot(x,ue,'-',x,u_rk,'.-',x,u_rkl,':',x,u_prk,'--',x,u_prrk,'.:');
hold on;
axis([-1,1,-0.1,1.1]);
set(gca,'xTick',[-1:0.4:1]);
set(gca,'yTick',[-0.1:0.2:1.1]);
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$u$','interpreter','latex','FontSize',16);
set(groot,'defaultLegendInterpreter','latex');
legend('Exact','RK(5,4),$\ \Delta t=0.002$','RK(5,4),\ $\Delta t=0.004385$','pRK(5,4),\ $\Delta t=0.5$','pRRK(5,4),\ $\Delta\hat{t}\approx0.009986$','location','NorthEastOutside');
set(gcf,'Units','centimeters','Position',[10 5 20 11],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);
