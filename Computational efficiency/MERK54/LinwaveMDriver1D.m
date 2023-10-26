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

% Solve Problem
time = 0; tstep = 0;
while (time<FinalTime)
  u = LinwaveM1D(x,u,h,k,kappa);
  load 'uex.mat';
  time = time+k;
  tstep = tstep+1;
end

toc
fprintf('Error_Linf=%e\n' ,norm(u-uex,inf))

figure(1)
plot(x,u,'r-o',x,uex,'k-.');
legend('numerical','exact')
axis([0,1,-1,1]);
axis square;
set(gcf,'color','w');