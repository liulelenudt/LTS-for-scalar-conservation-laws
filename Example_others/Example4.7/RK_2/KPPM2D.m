function [u] = KPPM2D(x,y,u,hx,hy,k,kappa);
% function [u] = KPPM2D(x,y,u,hx,hy,CFL,FinalTime);
% Purpose  : Integrate 2D KPP equation until FinalTime 
% using a monotone scheme.

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
psi(1) = 1;
psi(2) = psi(1)*( A(1,1)+k*kappa*b(1,1) );
psi(3) = psi(1)*( A(2,1)+k*kappa*b(2,1) ) + ...
         psi(2)*( A(2,2)+k*kappa*b(2,2) );
psi(4) = psi(1)*( A(3,1)+k*kappa*b(3,1) ) + ...
         psi(2)*( A(3,2)+k*kappa*b(3,2) ) + ...
         psi(3)*( A(3,3)+k*kappa*b(3,3) );
psi(5) = psi(1)*( A(4,1)+k*kappa*b(4,1) ) + ...
         psi(2)*( A(4,2)+k*kappa*b(4,2) ) + ...
         psi(3)*( A(4,3)+k*kappa*b(4,3) ) + ...
         psi(4)*( A(4,4)+k*kappa*b(4,4) );
psi(6) = psi(1)*( A(5,1)+k*kappa*b(5,1) ) + ...
         psi(2)*( A(5,2)+k*kappa*b(5,2) ) + ...
         psi(3)*( A(5,3)+k*kappa*b(5,3) ) + ...
         psi(4)*( A(5,4)+k*kappa*b(5,4) ) + ...
         psi(5)*( A(5,5)+k*kappa*b(5,5) );

maxvel=1;
rhsu = KPPMrhs2D(x,y,u,hx,hy,k,maxvel);
u1 = 1/psi(2)*  psi(1)*( A(1,1)*u  + b(1,1)*k*( rhsu  + kappa*u ) );
rhsu1 = KPPMrhs2D(x,y,u1,hx,hy,k,maxvel);
u2 = 1/psi(3)*( psi(1)*( A(2,1)*u  + b(2,1)*k*( rhsu  + kappa*u) )+...
                psi(2)*( A(2,2)*u1 + b(2,2)*k*( rhsu1 + kappa*u1 ) ) );
rhsu2 = KPPMrhs2D(x,y,u2,hx,hy,k,maxvel);
u3 = 1/psi(4)*( psi(1)*( A(3,1)*u  + b(3,1)*k*( rhsu  + kappa*u ) )+...
                psi(2)*( A(3,2)*u1 + b(3,2)*k*( rhsu1 + kappa*u1 ) )+...
                psi(3)*( A(3,3)*u2 + b(3,3)*k*( rhsu2 + kappa*u2 ) ) );
rhsu3 = KPPMrhs2D(x,y,u3,hx,hy,k,maxvel);
u4 = 1/psi(5)*( psi(1)*( A(4,1)*u  + b(4,1)*k*( rhsu  + kappa*u ) )+...
                psi(2)*( A(4,2)*u1 + b(4,2)*k*( rhsu1 + kappa*u1 ) )+...
                psi(3)*( A(4,3)*u2 + b(4,3)*k*( rhsu2 + kappa*u2 ) )+...
                psi(4)*( A(4,4)*u3 + b(4,4)*k*( rhsu3 + kappa*u3 ) ) );
rhsu4 = KPPMrhs2D(x,y,u4,hx,hy,k,maxvel);
u  = 1/psi(6)*( psi(1)*( A(5,1)*u  + b(5,1)*k*( rhsu  + kappa*u ) )+...
                psi(2)*( A(5,2)*u1 + b(5,2)*k*( rhsu1 + kappa*u1 ) )+...
                psi(3)*( A(5,3)*u2 + b(5,3)*k*( rhsu2 + kappa*u2 ) )+...
                psi(4)*( A(5,4)*u3 + b(5,4)*k*( rhsu3 + kappa*u3 ) )+...
                psi(5)*( A(5,5)*u4 + b(5,5)*k*( rhsu4 + kappa*u4 ) ) );
return