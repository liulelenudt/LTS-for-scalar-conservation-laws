function [u] = LinwaveM1D(x,u,h,k,kappa)
% MERK(5,4)
pz = @(z) 1+z+z^2/2+z^3/6+z^4/24+0.004477718303076*z^5;

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

ci = [0 0.39175222700392 0.58607968896779 0.47454236302687 0.93501063100924 1.0];
psi(1) = pz(ci(1)*k*kappa);
psi(2) = pz(ci(2)*k*kappa);
psi(3) = pz(ci(3)*k*kappa);
psi(4) = pz(ci(4)*k*kappa);
psi(5) = pz(ci(5)*k*kappa);
psi(6) = pz(ci(6)*k*kappa);
     
  % MERK(5,4)
  rhsu = LinwaveMrhs1D(x,u,h,k,1);
  u1 = 1/psi(2)*  psi(1)*( A(1,1)*u  + b(1,1)*k*( rhsu  + kappa*u ) );
  rhsu1 = LinwaveMrhs1D(x,u1,h,k,1);
  u2 = 1/psi(3)*( psi(1)*( A(2,1)*u  + b(2,1)*k*( rhsu  + kappa*u) )+...
                  psi(2)*( A(2,2)*u1 + b(2,2)*k*( rhsu1 + kappa*u1 ) ) );
  rhsu2 = LinwaveMrhs1D(x,u2,h,k,1);
  u3 = 1/psi(4)*( psi(1)*( A(3,1)*u  + b(3,1)*k*( rhsu  + kappa*u ) )+...
                  psi(2)*( A(3,2)*u1 + b(3,2)*k*( rhsu1 + kappa*u1 ) )+...
                  psi(3)*( A(3,3)*u2 + b(3,3)*k*( rhsu2 + kappa*u2 ) ) );
  rhsu3 = LinwaveMrhs1D(x,u3,h,k,1);
  u4 = 1/psi(5)*( psi(1)*( A(4,1)*u  + b(4,1)*k*( rhsu  + kappa*u ) )+...
                  psi(2)*( A(4,2)*u1 + b(4,2)*k*( rhsu1 + kappa*u1 ) )+...
                  psi(3)*( A(4,3)*u2 + b(4,3)*k*( rhsu2 + kappa*u2 ) )+...
                  psi(4)*( A(4,4)*u3 + b(4,4)*k*( rhsu3 + kappa*u3 ) ) );
  rhsu4 = LinwaveMrhs1D(x,u4,h,k,1);
  u  = 1/psi(6)*( psi(1)*( A(5,1)*u  + b(5,1)*k*( rhsu  + kappa*u ) )+...
                  psi(2)*( A(5,2)*u1 + b(5,2)*k*( rhsu1 + kappa*u1 ) )+...
                  psi(3)*( A(5,3)*u2 + b(5,3)*k*( rhsu2 + kappa*u2 ) )+...
                  psi(4)*( A(5,4)*u3 + b(5,4)*k*( rhsu3 + kappa*u3 ) )+...
                  psi(5)*( A(5,5)*u4 + b(5,5)*k*( rhsu4 + kappa*u4 ) ) );
return