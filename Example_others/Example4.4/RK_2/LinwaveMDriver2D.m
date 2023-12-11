% Driver script for solving the 2D Linwave equations using a monotone scheme
clear all
close all
clc
% Set problem parameters
Lx = 1; Ly = 1; Nx = 500; Ny = 500; hx = Lx/Nx; hy = Ly/Ny; 
FinalTime = 1.0;

% Define domain and initial conditions
xv = [0:hx:Lx]; yv = [0:hy:Ly]; [x,y] = meshgrid(xv,yv);
u = sin(4*pi*(x+y/2));
u0 = sin(4*pi*(x+y/2));
% Solve Problem
[u] = LinwaveM2D(x,y,u,hx,hy,FinalTime);

figure(1)
surf(xv,yv,u)
axis([0,1,0,1,-2,2]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[0:0.2:1]);
set(gca,'zTick',[-2:1.0:2]);
shading interp;
axis square;
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$y$','interpreter','latex','FontSize',16);
zlabel('$u$','interpreter','latex','FontSize',16);
set(gcf,'Units','centimeters','Position',[10 5 10 10],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
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
set(gcf,'Units','centimeters','Position',[10 5 12 10],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);

figure(3)
plot(xv,diag(u),'.-')
axis([0,1,-1,1]);
set(gca,'xTick',[0:0.2:1]);
set(gca,'yTick',[-1:0.5:1]);
shading interp;
axis square;
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$u$','interpreter','latex','FontSize',16);
set(gcf,'Units','centimeters','Position',[10 5 10 10],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);