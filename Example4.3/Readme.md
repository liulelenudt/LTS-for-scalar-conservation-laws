# Solving the 1D advection equation

# Description
Solving the 1D advection equation involves utilizing the SSPRK(5,4) scheme with a forward Euler time step, the SSPRK(5,4) scheme with a large time step, the pRK(5,4) scheme, and the pRRK(5,4) scheme for the SSP property.

# code
Driver code in LinwaveMDriver1D.m

Functions:

LinwaveLF - evaluate Lax-Friedrich numerical flux for the 1D advection equation.

LinwaveM1D - integrate the 1D advection equation until using a monotone scheme.

LinwaveMrhs1D - evaluate right hand side for the 1D advection equation using a monotone method.

relax - recalculate the last ralaxation time step using the pRRK(5,4) scheme.

wavetest - initial conditions for the 1D advection equation.

waveteste - exact solutions for the 1D advection equation.

extend - boundary conditions for the 1D advection equation.

# Example results
![image](https://github.com/liulelenudt/LTS-for-scalar-conservation-laws/assets/148626828/53bbbd3a-5e99-45a6-86d3-cd8630ecedc5)

Corresponding to Figure 3 in the manuscript:

(1) RK(5,4), $\Delta t=0.002$ corresponds to the SSPRK(5,4) scheme with a forward Euler time step, and $\Delta t=0.002$, $\kappa=0$;

(2) RK(5,4), $\Delta t=0.004385$ corresponds to the SSPRK(5,4) scheme with a large time step, and $\Delta t=0.004385$, $\kappa=0$;

(3) pRK(5,4), $\Delta t=0.5$ corresponds to the pRK(5,4) scheme, and $\Delta t=0.5$, $\kappa\approx496.984$;

(4) pRRK(5,4), $\Delta\hat{t}\approx0.009986$ corresponds to the pRRK(5,4) scheme, and $\Delta t=0.5$, $\Delta\hat{t}\approx0.009986$, $\kappa\approx496.984$.
