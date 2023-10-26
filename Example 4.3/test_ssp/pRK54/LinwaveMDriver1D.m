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

% pRK(5,4)
kappa = 1/h-1.508/k;
% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  [u] = LinwaveM1D(x,u,h,k,kappa);
  time = time+k;
  tstep = tstep+1;
end
toc

figure(1)
plot(x,u,'r-o',x,u0,'--');
legend('numerical','Exact')
axis([-1,1,-0.1,1.1]);
axis square;
set(gcf,'color','w');
