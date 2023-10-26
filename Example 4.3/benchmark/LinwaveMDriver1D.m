% Driver script for solving the 1D wave equation using a monotone scheme
clear all
close all
clc

% Set problem parameters
L = 2; 
FinalTime = 2; 
Nvec = 1000;
CFL = 0.90;
h = L/Nvec; 
k = 0.5;
% Define domain and initial conditions
x = [-1:h:1]'; 
[u] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
% [u0] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
[xe,uem] = extend(x,u,h,1,'P',0,'P',0);
u0 = uem(2:Nvec+2);

% pRK(3,3)
% c_s = (1.0+k*kappa+1.0/2.0*(k*kappa)^2)/...
%       (1.0+k*kappa+1.0/2.0*(k*kappa)^2+1.0/6.0*(k*kappa)^3);
% k_hat = c_s*k;

% pRK(5,4)
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
%   different "k" corresponding to fix(FinalTime/k_hat) "±1"
  [u] = LinwaveM1D(x,u,h,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u] = LinwaveM1D(x,u,h,k_res,kappa_res);
  break
  end
end

figure(1)
plot(x,u,'r-o',x,u0,'--');
legend('numerical','Exact')
axis([-1,1,-0.1,1.1]);
axis square;
set(gcf,'color','w');

figure(2)
load 'u_rk.mat';load 'u_rkl.mat';load 'u_prk.mat'
plot(x,u0,'-',x,u_rk,'.-',x,u_rkl,':',x,u_prk,'--',x,u,'.:');
hold on;
% plot(x,u,'o:','markersize',1.4)
axis([-1,1,-0.1,1.1]);
set(gca,'xTick',[-1:0.4:1]);
set(gca,'yTick',[-0.1:0.2:1.1]);
% axis square;
% set(gca,'FontName','Times New Roman','FontSize',12);
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$u$','interpreter','latex','FontSize',16);
% grid on;
set(groot,'defaultLegendInterpreter','latex');
legend('Exact','RK(5,4),$\ \Delta t=0.002$','RK(5,4),\ $\Delta t=0.004385$','pRK(5,4),\ $\Delta t=0.5$','pRRK(5,4),\ $\Delta\hat{t}\approx0.009986$','location','NorthEastOutside');
% 设置图片大小为13cm×10cm
set(gcf,'Units','centimeters','Position',[10 5 20 11],'color','w');
% 坐标线粗1.0磅
set(gca,'linewidth',1);
% Controls the box around the plotting area
set(gca,'box','on');
% 设置图中线宽1.2磅
set(get(gca,'Children'),'linewidth',1.2);
