% Driver script for solving the 1D wave equation using a monotone scheme
tic
clear
close all
clc

% Set problem parameters
L = 2; 
FinalTime = 2; 
Nvec = 1000;
h = L/Nvec; 
k = 0.5;
% Define domain and initial conditions
x = [-1:h:1]'; 
[u] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
[xe,uem] = extend(x,u,h,1,'P',0,'P',0);
u0 = uem(2:Nvec+2);

% pRRK(5,4)
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
%   different "k" corresponding to fix(FinalTime/k_hat) "¡À1"
  [u] = LinwaveM1D(x,u,h,k,kappa);
  time = time+k_hat;
  tstep = tstep+1;
  else
  [u] = LinwaveM1D(x,u,h,k_res,kappa_res);
  break
  end
end
toc

figure(1)
plot(x,u,'r-o',x,u0,'--');
legend('numerical','Exact')
axis([-1,1,-0.1,1.1]);
axis square;
set(gcf,'color','w');
