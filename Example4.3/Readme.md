# Solving the 1D advection equation

# Description
Solving the 1D advection equation involves utilizing the SSPRK(5,4) scheme with a forward Euler time step, the SSPRK(5,4) scheme with a large time step, the pRK(5,4) scheme, and the pRRK(5,4) scheme for the SSP problem.

The plots of SSPRK(5,4) scheme with a forward Euler time step, SSPRK(5,4) scheme with a large time step, pRK(5,4) scheme and pRRK(5,4) scheme correspond to Figure 3 in the manuscript

# code
Driver code in LinwaveMDriver1D.m

Functions:

LinwaveLF - evaluate Lax-Friedrich numerical flux for the 1D advection equation.

LinwaveM1D - integrate 1D wave equation until using a monotone scheme.

LinwaveMrhs1D - evaluate right hand side for the 1D advection equation using a monotone method.

relax - evaluate Lax-Friedrich numerical flux for the 1D advection equation.

wavetest - initial conditions for the 1D advection test problems.

waveteste - exact solution for the 1D advection test problems.

extend - boundary conditions for the 1D advection test problems.

# Example results

