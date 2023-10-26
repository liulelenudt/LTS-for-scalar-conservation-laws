% Driver script for solving the 2D Burgers equations using a monotone scheme
clear all
close all
clc
% Set problem parameters
Lx = 1; Ly = 1; Nx = 500; Ny = 500; hx = Lx/Nx; hy = Ly/Ny; 
FinalTime = 0.1; CFL = 0.9;
k = 0.01;%0.002835;
% Define domain and initial conditions
xv = [0:hx:Lx]; yv = [0:hy:Ly]; [x,y] = meshgrid(xv,yv);
u = sin(4*pi*(x+y/2));

% pRK(5,4)
kappa = (1/(hx/2)-1.508/k);
c_s = 1; % prk(5,4)
% c_s = (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+5*64/14293*(k*kappa)^4)/...
%       (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+1/24*(k*kappa)^4+64/14293*(k*kappa)^5);
k_hat = c_s*k;
% recalculate timestep of relaxation
[k_res,kappa_res] = relax(hx,k_hat,FinalTime);
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  if (time<=fix(FinalTime/k_hat)*k_hat)
%   different "k" corresponding to fix(FinalTime/k_hat) "±1"
  [u] = BurgersM2D(x,y,u,hx,hy,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u] = BurgersM2D(x,y,u,hx,hy,0.0013658599684643991165455432716917,0);
  break
  end
end

figure(1)
surf(xv,yv,u)
axis([0,1,0,1,-1,1]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[0:0.2:1]);
set(gca,'zTick',[-1:0.5:1]);
shading interp;
axis square;
% colorbar('ylim',[-1,1],'ytick',[-1:0.5:1])
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$y$','interpreter','latex','FontSize',16);
zlabel('$u$','interpreter','latex','FontSize',16);
% grid on;
% legend('Reference','RK(5,4),\Deltat=0.002','RK(5,4),\Deltat=0.004385','pRK(5,4),\Deltat=0.1','pRRK(5,4),\Deltat=0.1','location','West');
% 设置图片大小为12cm×10cm
set(gcf,'Units','centimeters','Position',[10 5 10 10],'color','w');
% 坐标线粗1.0磅
set(gca,'linewidth',1);
% Controls the box around the plotting area
set(gca,'box','on');
% 设置图中线宽1.2磅
set(get(gca,'Children'),'linewidth',1.2);


figure(2)
contourf(xv,yv,u,'LineStyle','none')
axis([0,1,0,1]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[0:0.2:1]);
shading interp;
axis square;
colorbar('ylim',[-1,1],'ytick',[-1:0.5:1])
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$y$','interpreter','latex','FontSize',16);
% grid on;
% legend('Reference','RK(5,4),\Deltat=0.002','RK(5,4),\Deltat=0.004385','pRK(5,4),\Deltat=0.1','pRRK(5,4),\Deltat=0.1','location','West');
% 设置图片大小为12cm×10cm
set(gcf,'Units','centimeters','Position',[10 5 12 10],'color','w');
% 坐标线粗1.0磅
set(gca,'linewidth',1);
% Controls the box around the plotting area
set(gca,'box','on');
% 设置图中线宽1.2磅
set(get(gca,'Children'),'linewidth',1.2);

figure(3)
plot(xv,diag(u),'.-')
axis([0,1,-1,1]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[-1:0.5:1]);
shading interp;
axis square;
% colorbar('ylim',[-1,1],'ytick',[-1:0.5:1])
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$u$','interpreter','latex','FontSize',16);
% grid on;
% legend('Reference','RK(5,4),\Deltat=0.002','RK(5,4),\Deltat=0.004385','pRK(5,4),\Deltat=0.1','pRRK(5,4),\Deltat=0.1','location','West');
% 设置图片大小为12cm×10cm
set(gcf,'Units','centimeters','Position',[10 5 10 10],'color','w');
% 坐标线粗1.0磅
set(gca,'linewidth',1);
% Controls the box around the plotting area
set(gca,'box','on');
% 设置图中线宽1.2磅
set(get(gca,'Children'),'linewidth',1.2);