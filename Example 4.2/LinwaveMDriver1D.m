% Driver script for solving the 1D wave equation using a monotone scheme
clear all
close all
clc

% Set problem parameters
L = 1; 
FinalTime = 0.5; 
Nvec = 2^8;
dt = 2^(-9);
kvec = [dt*8 dt*4 dt*2 dt dt/2 dt/4 dt/8 dt/16 dt/32 dt/64 dt/128 dt/1024];
CFL = 0.90;
Error0 = []; Error1 = []; Error2 = [];
order0 = []; order1 = []; order2 = [];

for i = 1:length(kvec)
    k = kvec(i);
    h = L/Nvec; 
    
    % Define domain and initial conditions
    x = [0:h:1]'; 
    % [u] = wavetest(x,0.5,-0.7,0.005,10,log(2)/(36*0.005^2));
    [u] = sin(2*pi*x);
    [uex] = sin(2*pi*x);

    % Solve Problem
    [u] = LinwaveM1D(x,u,h,k,FinalTime);
    [uex] = LinwaveM1D(x,uex,h,2^(-18),FinalTime);
    error0 = norm(u-uex,inf);
    error1 = norm(u-uex,1)*h;
    error2 = norm(u-uex,2)*sqrt(h);
    Error0 = [Error0,error0];
    Error1 = [Error1,error1];
    Error2 = [Error2,error2];
end
% computating convergence order
for n = 1:length(kvec)-1
    order0(n+1) = log(Error0(n)/Error0(n+1))/(log(kvec(n)/kvec(n+1)));
    order1(n+1) = log(Error1(n)/Error1(n+1))/(log(kvec(n)/kvec(n+1)));
    order2(n+1) = log(Error2(n)/Error2(n+1))/(log(kvec(n)/kvec(n+1)));
end
fprintf('***************************************************************\n')
fprintf(' Accuracy of temporal discritization\n');
fprintf('***************************************************************\n')
fprintf(' nT \t Linf-Norm \t Degree \t L2-Norm \t Degree\n');
for n = 1:length(kvec)
    fprintf('%3.0f \t %1.2e \t %2.2f \t\t %1.2e \t %2.2f \n',...
        1/kvec(n)/2^9,Error0(n),order0(n),Error2(n),order2(n));
end

figure(1)
plot(x,u,'r-o',x,uex,'k-.');
legend('numerical','exact')
axis([0,1,-1,1]);
axis square;
set(gcf,'color','w');

figure(2)
load 'Err_prrk.mat'; load 'Err_prk.mat'; load 'Err_ifrk.mat';
loglog(kvec,Error0(:),'-s',kvec,Err_ifrk(:),'->',kvec,Err_prk(:),'-<',kvec,Err_prrk(:),'-o')
hold on;
loglog(kvec, 15000000*kvec.^4, 'k--', 'linewidth', 1);
annotation('textarrow', [0.46 0.42], [0.31, 0.32], 'string', '$\Delta t^4$', 'interpreter','latex','fontsize', 12);
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$\Delta t$','interpreter','latex','FontSize',16);
ylabel('$L^{\infty}~Error$','interpreter','latex','FontSize',16);
set(gcf,'Units','centimeters','Position',[10 1 15 15],'color','w');
legend('RK(5, 4)','IFRK(5, 4)','pRK(5, 4)','pRRK(5, 4)','location','NorthWest');
axis([0.5e-5,kvec(1),0.5e-13,2.5e-0]);
set(gca,'linewidth',1);
% Controls the box around the plotting area
set(gca,'box','on');
set(get(gca,'Children'),'linewidth',1.2);

