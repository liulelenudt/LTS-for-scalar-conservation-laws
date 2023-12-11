function [k_res,kappa_res] = relax(h,k_hat,FinalTime)
syms k_res
kappa_res = 1/h-1.508/k_res;
k_res = vpasolve((1+k_res*kappa_res+1/2*(k_res*kappa_res)^2+1/6*(k_res*kappa_res)^3+5*64/14293.0*(k_res*kappa_res)^4)/...
                 (1+k_res*kappa_res+1/2*(k_res*kappa_res)^2+1/6*(k_res*kappa_res)^3+1/24*(k_res*kappa_res)^4+64/14293*(k_res*kappa_res)^5)-...
                 (FinalTime-fix(FinalTime/k_hat)*k_hat)/k_res);
k_res = k_res(1);
kappa_res = max(1/h-1.508/k_res,0);
end