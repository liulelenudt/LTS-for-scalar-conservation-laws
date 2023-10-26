% u'(t)=1-c*sqrt(u)*u, t in [0,0.05] 
% Initial condition: u(0)=u_s ;
clear; clc,close all;
tic
Nt = [20 40 80];             % Number of partitions
U = [];
Error = [];
c = 10000;
u_s = 1/sqrt(c);
ta = 0;
tb = 0.05;
kappa = 800;
fun = @(t,u) 1-c*abs(u)*u;                % RHS
pz = @(z) 1+z+z^2/2+z^3/6+z^4/24+0.004477718303076*z^5;
for k = 1:length(Nt)
    N = Nt(k);
    h = (tb-ta)/N;    
    % relaxation
    c_s = 1;
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
    u(1) = 0.9*u_s;                       % initial value
    
    ci = [0 0.39175222700392 0.58607968896779 0.47454236302687 0.93501063100924 1.0];
    psi(1) = pz(ci(1)*h*kappa);
    psi(2) = pz(ci(2)*h*kappa);
    psi(3) = pz(ci(3)*h*kappa);
    psi(4) = pz(ci(4)*h*kappa);
    psi(5) = pz(ci(5)*h*kappa);
    psi(6) = pz(ci(6)*h*kappa);
    for n = 1:N_r
        % MERK(5,4)
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
figure(1)
plot(T{1},U{1},'-',T{2},U{2},':.',T{3},U{3},'--',t_s,u_s,':');
axis([ta,tb,0.0085,0.0110]);
set(gca,'xTick',[0:0.01:0.05]);
set(gca,'yTick',[0.0085:0.0005:0.0110]);
% get(gca,'xtick') ; % 得到坐标的实际大小
set(gca,'xticklabel',get(gca,'xtick')); % 将x显示的字符替换为实际大小
% get(gca,'ytick');
set(gca,'yticklabel',get(gca,'ytick'));

set(gcf,'color','w');
set(gca,'FontSize',12,'fontweight','demi');
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('$u$','Interpreter','latex','FontSize',16);
% title('$ $','Interpreter','latex','FontSize',16,'FontWeight','bold');
legend('$N=20$','$N=40$','$N=80$','Equilibrium','location','southeast');
set(gcf,'Units','centimeters','Position',[10 5 12 9]);%设置图片大小为13cm×10cm
set(gca,'linewidth',1); %坐标线粗0.5磅
set(gca,'box','off');%Controls the box around the plotting area
set(get(gca,'Children'),'linewidth',2);%设置图中线宽1磅
% box on; %封闭边框

figure(2)
semilogy(T{1},abs(U{1}-0.01),'-',T{2},abs(U{2}-0.01),':.',T{3},abs(U{3}-0.01),'--',t_s,u_s-u_s,':');
axis([ta,tb,1e-8,1e-0]);
% get(gca,'xtick') ; % 得到坐标的实际大小
% set(gca,'xticklabel',get(gca,'xtick')); % 将x显示的字符替换为实际大小
% get(gca,'ytick');
% set(gca,'yticklabel',get(gca,'ytick'));

set(gcf,'color','w');
set (gca,'FontSize',12,'fontweight','demi');
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('$|u-u^{*}|$','Interpreter','latex','FontSize',16);
% title('$ $','Interpreter','latex','FontSize',16,'FontWeight','bold');
legend('$N=20$','$N=40$','$N=80$','Equilibrium','location','northeast');
set(gcf,'Units','centimeters','Position',[10 5 13 10]);%设置图片大小为13cm×10cm
set(gca,'linewidth',1); %坐标线粗0.5磅
set(gca,'box','off');%Controls the box around the plotting area
set(get(gca,'Children'),'linewidth',2);%设置图中线宽1磅
toc
