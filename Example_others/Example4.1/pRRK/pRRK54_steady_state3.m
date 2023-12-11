% u'(t)=1-c*sqrt(u)*u, t in [0,0.05] 
% Initial condition: u(0)=u_s ;
clear; clc,close all;
tic
Nt = [20 40 80];
U = [];
Error = [];
c = 10000;
u_s = 1/sqrt(c);
ta = 0;
tb = 0.05;
kappa = 800;
fun = @(t,u) 1-c*abs(u)*u;
for k = 1:length(Nt)
    N = Nt(k);
    h = (tb-ta)/N;    
    % relaxation
    c_s = (1+h*kappa+1/2*(h*kappa)^2+1/6*(h*kappa)^3+5*64/14293*(h*kappa)^4)/...
          (1+h*kappa+1/2*(h*kappa)^2+1/6*(h*kappa)^3+1/24*(h*kappa)^4+64/14293*(h*kappa)^5);
    h_r = c_s*h;
    N_r = (tb-ta)/h_r;
    A = [1 0 0 0 0;
         0.444370493651235 0.555629506348765 0 0 0;
         0.620101851488403 0 0.379898148511597 0 0;
         0.178079954393132 0 0 0.821920045606868 0;
         0 0 0.517231671970585 0.096059710526147 0.386708617503268];
    b = [0.39175222657189 0 0 0 0;
         0 0.368410593050371 0 0 0;
         0 0 0.251891774271694 0 0;
         0 0 0 0.54497475022852 0;
         0 0 0 0.063692468666290 0.226007483236906];
    t = ta:h_r:tb;                   % interval partition
    u(1) = 1.1*u_s;                  % initial value
    % recursive approximations
    psi(1) = 1;
    psi(2) = psi(1)*( A(1,1)+h*kappa*b(1,1) );
    psi(3) = psi(1)*( A(2,1)+h*kappa*b(2,1) ) + ...
             psi(2)*( A(2,2)+h*kappa*b(2,2) );
    psi(4) = psi(1)*( A(3,1)+h*kappa*b(3,1) ) + ...
             psi(2)*( A(3,2)+h*kappa*b(3,2) ) + ...
             psi(3)*( A(3,3)+h*kappa*b(3,3) );
    psi(5) = psi(1)*( A(4,1)+h*kappa*b(4,1) ) + ...
             psi(2)*( A(4,2)+h*kappa*b(4,2) ) + ...
             psi(3)*( A(4,3)+h*kappa*b(4,3) ) + ...
             psi(4)*( A(4,4)+h*kappa*b(4,4) );
    psi(6) = psi(1)*( A(5,1)+h*kappa*b(5,1) ) + ...
             psi(2)*( A(5,2)+h*kappa*b(5,2) ) + ...
             psi(3)*( A(5,3)+h*kappa*b(5,3) ) + ...
             psi(4)*( A(5,4)+h*kappa*b(5,4) ) + ...
             psi(5)*( A(5,5)+h*kappa*b(5,5) );
         
    for n = 1:N_r
        % pRRK(5,4) scheme
        u1 = 1/psi(2)*( psi(1)*( A(1,1)*u(n)+b(1,1)*h.*( fun(t(n),u(n))+       kappa*u(n) ) ) );
        u2 = 1/psi(3)*( psi(1)*( A(2,1)*u(n)+b(2,1)*h.*( fun(t(n),u(n))+       kappa*u(n) ) )+...
                        psi(2)*( A(2,2)*u1+  b(2,2)*h.*( fun(t(n)+0.39175222700392*h_r,u1)+kappa*u1 ) ) );
        u3 = 1/psi(4)*( psi(1)*( A(3,1)*u(n)+b(3,1)*h.*( fun(t(n),u(n))+       kappa*u(n) ) )+...
                        psi(2)*( A(3,2)*u1+  b(3,2)*h.*( fun(t(n)+0.39175222700392*h_r,u1)+kappa*u1 ) )+...
                        psi(3)*( A(3,3)*u2+  b(3,3)*h.*( fun(t(n)+0.58607968896779*h_r,u2)+kappa*u2 ) ) );
        u4 = 1/psi(5)*( psi(1)*( A(4,1)*u(n)+b(4,1)*h.*( fun(t(n),u(n))+       kappa*u(n) ) )+...
                        psi(2)*( A(4,2)*u1+  b(4,2)*h.*( fun(t(n)+0.39175222700392*h_r,u1)+kappa*u1 ) )+...
                        psi(3)*( A(4,3)*u2+  b(4,3)*h.*( fun(t(n)+0.58607968896779*h_r,u2)+kappa*u2 ) )+...
                        psi(4)*( A(4,4)*u3+  b(4,4)*h.*( fun(t(n)+0.47454236302687*h_r,u3)+kappa*u3 ) ) );
    u(n+1) = 1/psi(6)*( psi(1)*( A(5,1)*u(n)+b(5,1)*h.*( fun(t(n),u(n))+       kappa*u(n) ) )+...
                        psi(2)*( A(5,2)*u1+  b(5,2)*h.*( fun(t(n)+0.39175222700392*h_r,u1)+kappa*u1 ) )+...
                        psi(3)*( A(5,3)*u2+  b(5,3)*h.*( fun(t(n)+0.58607968896779*h_r,u2)+kappa*u2 ) )+...
                        psi(4)*( A(5,4)*u3+  b(5,4)*h.*( fun(t(n)+0.47454236302687*h_r,u3)+kappa*u3 ) )+...
                        psi(5)*( A(5,5)*u4+  b(5,5)*h.*( fun(t(n)+0.93501063100924*h_r,u4)+kappa*u4 ) ) );
    end
    T{k} = t;
    U{k} = u;
end
t_s = ta:0.0001:tb;
u_s = ones(1,length(t_s))*u_s;

%% logarithmic plot
figure(1)
semilogy(T{1},abs(U{1}-0.01),'-',T{2},abs(U{2}-0.01),':.',T{3},abs(U{3}-0.01),'--');
axis([ta,tb,1e-8,0.01]);
set(gcf,'color','w');
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('$|u-u^{*}|$','Interpreter','latex','FontSize',16);
legend('$N=20$','$N=40$','$N=80$','location','northeast');
set(gcf,'Units','centimeters','Position',[10 5 13 10]);
set(gca,'linewidth',1);
set(gca,'box','off');
set(get(gca,'Children'),'linewidth',2);
toc