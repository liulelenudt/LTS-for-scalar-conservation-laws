% Driver script for solving the 2D Linwave equations using a monotone scheme
clear all
close all
clc
% Set problem parameters
Lx = 1; Ly = 1; Nx = 500; Ny = 500; hx = Lx/Nx; hy = Ly/Ny; 
FinalTime = 1.0; CFL = 0.9;

% Define domain and initial conditions
xv = [0:hx:Lx]; yv = [0:hy:Ly]; [x,y] = meshgrid(xv,yv);
u = sin(4*pi*(x+y/2));
u0 = sin(4*pi*(x+y/2));
% Solve Problem
[u] = LinwaveM2D(x,y,u,hx,hy,CFL,FinalTime);

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