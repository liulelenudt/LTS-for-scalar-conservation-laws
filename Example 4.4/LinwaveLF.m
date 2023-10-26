function [numflux] = LinwaveLF(u,v,lambda,maxvel)
% function [numflux] = LinwaveLF(u,v,lambda,maxvel);
% Purpose: Evaluate the Lax Friedrich numerical flux for Linwave equation

fu = u; fv = v;
numflux = (fu+fv)/2 - maxvel/2*(v-u);
end