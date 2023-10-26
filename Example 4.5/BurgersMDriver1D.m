% Driver script for solving the 1D Burgers equation using a monotone scheme
clear all
close all
clc
% Set problem parameters
L = 1; FinalTime = 0.4; N = 500; h = L/N; CFL = 0.9;

% Define domain and initial conditions
x = [0:h:1]'; 
k = 0.01;
u = sin(2*pi*x); %periodic BC needed
u0 = sin(2*pi*x);
% [u] = wavetest(x);
% [u0] = wavetest(x);

% u = (1-sign(x-0.2))/2+1; % Constant BC needed

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
  [u] = BurgersM1D(x,u,h,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u] = BurgersM1D(x,u,h,k_res,kappa_res);
  break
  end
end

for i = 1:N+1
    x(i) = 0 + (i-1)*h;
%     xn = fixpoint(x(i),t);
    xn = bisection(x(i),FinalTime);
    ue(i) = sin(2*pi*xn);
end

figure(1)
plot(x,u,'-o',x,ue,'-s');
legend('numerical','exact')
% axis([0,1,0,1]);
axis square;
set(gcf,'color','w');

figure(2)
% load 'u_rk1.mat';load 'u_rkl1.mat';load 'u_prk1.mat';
% plot(x,ue,'-',x,u_rk1,'.-',x,u_rkl1,':',x,u_prk1,'--',x,u,'.:');
load 'u_rk2.mat';load 'u_rkl2.mat';load 'u_prk2.mat';
plot(x,ue,'-',x,u_rk2,'.-',x,u_rkl2,':',x,u_prk2,'--',x,u,'.:');
hold on;
axis([0,1,-1.2,1.2]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[-1.2:0.4:1.2]);
axis square;
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$u$','interpreter','latex','FontSize',16);
% grid on;
legend('Exact','RK(5,4),\ $\Delta t=0.002$','RK(5,4),\ $\Delta t=0.005395$','pRK(5,4),\ $\Delta t=0.01$','pRRK(5,4),\ $\Delta\hat{t}\approx0.008019$','location','Northeast');
% 设置图片大小为13cm×10cm
set(gcf,'Units','centimeters','Position',[10 1 15.6 15],'color','w');
% 坐标线粗1.0磅
set(gca,'linewidth',1);
% Controls the box around the plotting area
set(gca,'box','on');
% 设置图中线宽1.2磅
set(get(gca,'Children'),'linewidth',1.2);
