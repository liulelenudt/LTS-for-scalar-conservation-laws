% Driver script for solving the 2D KPP equations using a monotone scheme
clear all
close all
clc
% Set problem parameters
Nx = 500; Ny = 500; hx = 4.0/Nx; hy = 4.0/Ny; 
FinalTime = 1.0;

% Define domain and initial conditions
xv = [-2:hx:2]; yv = [-2.5:hy:1.5]; [x,y] = meshgrid(xv,yv);
u = pi/4 + (sqrt(x.^2+y.^2)<=1)*13*pi/4;

%% RK(5,4) scheme with large time step
k = 0.009895;
kappa = 0;

% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u] = KPPM2D(x,y,u,hx,hy,k,kappa);
  time = time+k;
  tstep = tstep+1;
end

%% plots
figure(1)
surf(xv,yv,u)
axis([-2,2,-2.5,1.5,0,11]);
set(gca,'xTick',[-2:1:2]);
set(gca,'yTick',[-2.5:1:1.5]);
set(gca,'zTick',[0:2:10]);
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
axis([-2,2,-2.5,1.5]);
set(gca,'xTick',[-2:1:2]);
set(gca,'yTick',[-2.5:1:1.5]);
shading interp
axis square;
colorbar('ylim',[0,11],'ytick',[0:2:11])
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$y$','interpreter','latex','FontSize',16);
set(gcf,'Units','centimeters','Position',[10 5 12 10],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);

figure(3)
plot(xv,flip(diag(fliplr(u))),'.-')
axis([-2,2,0,11.5]);
set(gca,'xTick',[-2:1:2]);
set(gca,'yTick',[0:2:11.5]);
shading interp;
axis square;
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$x$','interpreter','latex','FontSize',16);
ylabel('$u$','interpreter','latex','FontSize',16);
set(gcf,'Units','centimeters','Position',[10 5 10 10],'color','w');
set(gca,'linewidth',1);
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);
