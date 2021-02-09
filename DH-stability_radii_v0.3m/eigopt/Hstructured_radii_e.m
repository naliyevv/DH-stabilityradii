function [f, g, t, tsdp, tchol] = Hstructured_radii_e(omega,pars)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% Auxiliary routine called by structuredradiiJR_str, HstructuredradiiJR_ss, structuredradiiJR_Qinv_str
%
% TASK:
% This computes the objective eigenvalue function (f) and its derivative (g) at a given omega
% to be minimized for the estimation of the structured stability radius.
%
% The quantity computed (f) is tilde{eta}^{Herm}(R; B, i*omega) defined in eqn (5.1) in [1].
% This corresponds to the square of the backward error of i*omega to be an eigenvalue of (J-R)Q
% subject to structured perturbations of R of the form B Delta B^*
%
% J, Q, B, (J-R)Q must be passed through pars.J, pars.Q, pars.B, pars.A.
%
% REFERENCES:
% [1] - Computation of Stability Radii for Large-Scale Port-Hamiltonian Systems. Tech. Report.
%       Nicat Aliyev, Volker Mehrmann and Emre Mengi.



n = size(pars.J,1);
m = size(pars.B,2);
h = 10^-8;
Q = pars.Q;
J = pars.J;
B = pars.B;









%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%% (1) Compute f (i.e, the backward error at omega)
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%


                                   
                                   
%%%%%%%%%%%%%%%%%
%%%% Compute Z = ((J - R)Q - i*omega I)^{-1} B
%%%% H0, H1 are defined in terms of Z.
%%%%%%%%%%%%%%%%%
tic;
if issparse(pars.J)
    W = pars.A - 1i*omega*speye(n);
else
    W = pars.A - 1i*omega*eye(n);
end
Z = W\B;



%%%%%%%%%%%%%
% FORM H0 and THE CHOLESKY FACTOR
%%%%%%%%%%%%%
H0 = B'*Q*Z;
H0 = H0'*H0;
H0 = (H0 + H0')/2;

% should normally not happen, but in case of rounding errors
if (min(real(eig(H0))) < 10^-14)
    H0 = H0 + (-min(real(eig(H0)))+10^-14)*eye(m);
end
      
% L = chol(H0,'lower');
 R = chol(H0);
 L = R';
      
H0 = L\(L'\eye(m));
H0 = (H0 + H0')/2;
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%


%%%%%%%%%%%%%
% FORM H1
%%%%%%%%%%%%%
H1 = L\(Z'*(Q*(B/L')));
H1 = 1i*(H1 - H1');
H1 = (H1 + H1')/2;
tchol = toc;
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
      
      
%%%%%%%%%%%%%%
% In case H1 is definite, the eigenvalue i*omege is not attainable
% with Hermitian perturbations of R of the form B Delta B^*.
% In theory f is infinity, but we set it to some large value.
%%%%%%%%%%%%%%
if (min(eig(H1)) > 0) || (max(eig(H1)) < 0)
      f = 0.1;
      g = 0;
      t = inf;
      tsdp = 0;
      return;
end
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%

      
%%%%%%%%%%%%%%
% Maximize lambda_min(H0 + t H1) over t using eigopt
% (Alternatively an SDP solver can be used for this purpose, but this seems much slower)
%%%%%%%%%%%%%%
parsin.H0 = H0;
parsin.H1 = H1;
parsin.tol = 10^-12;
parsin.gamma = -0.0000000001;
parsin.minmax = 1;
intvalin.lb = -20;
intvalin.ub = 20;

      
tic;
[f,t,parsout] = eigopt('Hstructured_radii_inner_e',intvalin,parsin);
tsdp = toc;
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
      
      
if (f >= 0.1)
      f = 0.1;
      g = 0;
      t = inf;
      return;
end

   
      
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%% End of (1)
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
      

      
      
      
      
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For now we use finite differences to estimate the derivative (g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

      
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%% (2) Now compute the backward error at (omega+h), this is for finite differences
%%%     Same calculations as in (1) but at (omega+h) instead of omega
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
      
      
tic;
if issparse(pars.J)
    W = pars.A - 1i*(omega+h)*speye(n);
else
    W = pars.A - 1i*(omega+h)*eye(n);
end
Z = W\B;


% FORM H0 and THE CHOLESKY FACTOR
H0 = B'*Q*Z;
H0 = H0'*H0;
H0 = (H0 + H0')/2;
R = chol(H0);
L = R';
      
H0 = L\(L'\eye(m));
H0 = (H0 + H0')/2;

% FORM H1
H1 = L\(Z'*(Q*(B/L')));
H1 = 1i*(H1 - H1');
H1 = (H1 + H1')/2;
tchol2 = toc;

parsin.H0 = H0;
parsin.H1 = H1;
      
tic;
[fh,th] = eigopt('Hstructured_radii_inner_e',intvalin,parsin);
tsdp2 = toc;

      
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%% End of (2)
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
      

      
      
g = (fh - f)/h;

      
tchol = tchol + tchol2;
tsdp = tsdp + tsdp2;

return;
