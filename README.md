# DH-stabilityradii
Matlab implementations of the algorithms tailored in  to compute the non-Hermitian and Hermitian restricted stability radii for a linear dissipative Hamiltonian (DH) system.
Usage:
The main routines to call are the following.


non-Hermitian Stability Radii 
DHradiiJR_nonHermit  : This is for the large-scale computation of the non-Hermitian stability radius with respect to perturbations of J and R.
DHradiiQ_nonHermit  : This is for the large-scale computation of the non-Hermitian stability radius with respect to perturbations of Q.
DHradiiJR_nonHermit_Qinv  : This is similar to DHradiiJR_nonHermit , but takes the inverse of Q as an input (instead of Q). Use this rather than DHradiiJR_nonHermit if the inverse of Q is available and sparse rather than Q (which may potentially be dense).

Hermitian Stability Radii
DHradiiJR_Hermit_ss : For the small-scale computation of the Hermitian stability radius with respect to Hermitian perturbations of R.
DHradiiJR_Hermit : For the large-scale computation of the Hermitian stability radius with respect to Hermitian perturbations of R.
DHradiiJR_Hermit_Qinv : Similar to DHradiiJR_Hermit , however expects the inverse of Q as an input rather than Q. Use this if the inverse of Q is available and sparse.

For the usage of these routines please see the sample calls. The input and output arguments for these routines are also described at the beginning of the routines.
