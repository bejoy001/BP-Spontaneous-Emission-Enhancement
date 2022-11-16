This repository contains the required MATLAB codes used to generate the figures in paper 
"Deep-Utraviolet-Spontaneous-Emission-Enhanced-by-Layer-Dependent-Black-Phosphorus-Plasmonics" which is under review in
"Optical Express". The doi link will be shared after complete acceptance.

The codes have been written in MATLAB. The folder contains three .m files:

1. BP_Optical_Conductivity.m computes the optical conductivity of Black Phosphorus using Kubo Formula given in
equation 10 in the paper.
2. Purcell_Factor.m calculates Purcell Factor using generalized dyadic Green's Function Formalism which require the
optical conductivity calculated using BP_Optical_Conductivity.m
3. Finally Rate_Equation.m solves the rate equations for different possible transitions of Hydrogen atom to compute
the spectrum. 

