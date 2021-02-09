function [f, z, parsout] = DHradiiJR_Hermit_ss(J,Q,R,B,intval)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% TASK:
% Computes the structured stability radius for the
% dissipative-Hamiltonian (DH) system of the form
%                       x' =(J-R)Qx     (1)
% with respect to Hermitian pertubations of R.
%
% This is meant only for small-scale problems, and does not use subspace projections.
%
% Recall that the structured stability radius is defined by
%
%    inf { || Delta || | Delta = Delta' and (J -  (R + B Delta B'))Q has an eigenvalue on the imaginary axis  }
%
% for a given full-rank tall-skinny matrix B.
%
% The approach exploits the eigenvalue optimization characterization in [1, Theorem 3.2].
%
%
% CALL : [f, z, parsout] = DHradiiJR_ss(J,Q,R,B,intval)
%
%
% INPUT:
%     J: skew-symmetric matrix for the energy flux of the system
%     Q: symmetric, positive definite matrix for the total energy of the system
%     R: symmetric positive semi-definite dissipation matrix of the system
%     B: the restriction matrix (must be full-rank)
%     intval: specifies the interval where initial interpolation is to be performed
%             with the fields
%               intval.l, intval.u
%               the interval for initial interpolation is [i*intval.l, i*intval.u]
%     (default - [intval.l, intval.u] = [-2, 0 ])
%
%
% OUTPUT:
%       f: computed structured stability radius of the DH system, i.e.,
%       z: the real number such that i*z becomes an eigenvalue under the smallest Hermitian perturbation
%          possible of norm equal to the structured stability radius.
%    parsout: structure containing information regarding the progress of the algorithm
%           parsout.time: total time consumed by the algorithm to reach prescribed accuracy
%
%
% REFERENCES:
% [1] - Computation of Stability Radii for Large-Scale Port-Hamiltonian Systems. Tech. Report.
%       Nicat Aliyev, Volker Mehrmann and Emre Mengi.




time1 = cputime;


if (nargin < 5)
    intval.l = -2;
    intval.u = 0;
end




%%%%%%%%%%%%
% Compute directly using eigopt
%%%%%%%%%%%%
pars.J = J;
pars.R = R;
pars.Q = Q;
pars.A = (J - R)*Q;
pars.B = B;


pars.tol = 10^-6;
pars.gamma = -1;

intval.lb = intval.l;
intval.ub = intval.u;
pars.isprint = 0;
pars.itertol = 100000;

[f,z,parsout] = eigopt('Hstructured_radii_e',intval,pars);

%%%%%%%%%%%%%
%%%%%%%%%%%%%


f = sqrt(f);

time2 = cputime;


parsout.time = time2 - time1;


return;
