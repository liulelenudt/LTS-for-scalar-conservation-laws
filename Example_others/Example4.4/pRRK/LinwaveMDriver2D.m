% Driver script for solving the 2D Linwave equations using a monotone scheme
clear all
close all
clc
% Set problem parameters
Lx = 1; Ly = 1; Nx = 500; Ny = 500; hx = Lx/Nx; hy = Ly/Ny; 
FinalTime = 1.0;
k = 0.1;
% Define domain and initial conditions
xv = [0:hx:Lx]; yv = [0:hy:Ly]; [x,y] = meshgrid(xv,yv);
[u] = sin(4*pi*(x+y/2));
u0 = sin(4*pi*(x+y/2));

% pRRK(5,4)
kappa = (1/(hx/2)-1.508/k);
c_s = (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+5*64/14293*(k*kappa)^4)/...
      (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+1/24*(k*kappa)^4+64/14293*(k*kappa)^5);
k_hat = c_s*k;
% recalculate timestep of relaxation
[k_res,kappa_res] = relax(hx,k_hat,FinalTime);
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  if (time<=fix(FinalTime/k_hat)*k_hat)
  [u] = LinwaveM2D(x,y,u,hx,hy,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u] = LinwaveM2D(x,y,u,hx,hy,0.0036405853698093303698550329355368,171.56179744853558929770433333916);
  % k_res=0.0036405853698093303698550329355368;  kappa_res=171.56179744853558929770433333916
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