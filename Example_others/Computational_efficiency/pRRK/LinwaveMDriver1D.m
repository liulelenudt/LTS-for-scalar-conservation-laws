% Driver script for solving the 1D wave equation using a monotone scheme
clear
close all
clc
tic
% Set problem parameters
L = 1; 
FinalTime = 2; 
Nvec = 2^12;
k = 2^(-10);
h = L/Nvec; 

% Define domain and initial conditions
x = [0:h:1]'; 
u = sin(2*pi*x);
% [uex] = sin(2*pi*x);

% Solve Problem
% kappa = 0;
kappa = 1/h;
% c_s = 1;
c_s = (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+5*64/14293*(k*kappa)^4)/...
      (1+k*kappa+1/2*(k*kappa)^2+1/6*(k*kappa)^3+1/24*(k*kappa)^4+64/14293*(k*kappa)^5);
k_hat = c_s*k;
% recalculate timestep of relaxation
[k_res,kappa_res] = relax(h,k_hat,FinalTime);
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  if (time<=fix(FinalTime/k_hat-1)*k_hat)
  u = LinwaveM1D(x,u,h,k,kappa);
  load 'uex.mat';
  time = time+k_hat;
  tstep = tstep+1;
  else
  u = LinwaveM1D(x,u,h,k_res,kappa_res);
  break
  end
end

toc
fprintf('Error_Linf=%e\n' ,norm(u-uex,inf))

figure(1)
plot(x,u,'r-o',x,uex,'k-.');
legend('numerical','exact')
axis([0,1,-1,1]);
axis square;
set(gcf,'color','w');