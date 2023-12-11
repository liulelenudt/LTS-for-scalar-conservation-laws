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

%% pRRK(5,4) scheme
k = 0.1;
kappa = (1/(hx/2)-1.508/k);
c_s = (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+5*64/14293*(k*kappa)^4)/...
      (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+1/24*(k*kappa)^4+64/14293*(k*kappa)^5);
k_hat = c_s*k;
% recalculate timestep of relaxation
[k_res,kappa_res] = relax(hx,k_hat,FinalTime);
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  if (time<=fix(FinalTime/k_hat-1)*k_hat)
  [u] = KPPM2D(x,y,u,hx,hy,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u] = KPPM2D(x,y,u,hx,hy,0.016741081175214379267407323695646,319.84437155317847215717449646672);
  % k_res=0.016741081175214379267407323695646;  kappa_res=319.84437155317847215717449646672
  break
  end
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
