function [f, g, t, tsdp, tchol] = Hstructured_radii_e_invQs(omega,pars,modeQ,modeQs,mylu)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% Auxiliary routine called by structuredradiiJR_Qinv_str
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
% ************************
% *** IMPORTANT REMARK ***
% ************************
% USE THIS IF Q^{-1} IS AVAILABLE AND SPARSE RATHER THAN Q.
% OTHERWISE IT IS BETTER TO USE Hstructured_radii_e.m
% ************************
% ************************
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
R = pars.R;
                                   








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

if modeQ

    if issparse(pars.J)
       W = pars.A - 1i*omega*speye(n);
    else
       W = pars.A - 1i*omega*eye(n);
    end
    Z = W\B;
                                   
else
                                   
    Qinv = pars.Qinv;
    LQinv = pars.LQinv;
    UQinv = pars.UQinv;
    PQinv = pars.PQinv;
                                   
    if ~modeQs
       Z = Qinv*(((J-R)-omega*1i*Qinv)\B);
    else
       pars.s = 1i*omega;
       pars.B = B;
       pars.Qinv = Qinv;
       pars.sq = 0;
       [Z,~] = feval(mylu,(J-R)-omega*1i*Qinv,pars);
                                   
    end
end
                                   
                                   

    
                                   
%%%%%%%%%%%%%
% FORM H0 and THE CHOLESKY FACTOR
%%%%%%%%%%%%%
H0 = B'* (UQinv \ (LQinv \ (PQinv*Z)));
H0 = H0'*H0;
H0 = (H0 + H0')/2;
% L = chol(H0,'lower');
 Rc = chol(H0);
 L = Rc';
      
H0 = L\(L'\eye(m));
H0 = (H0 + H0')/2;
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%


%%%%%%%%%%%%%
% FORM H1
%%%%%%%%%%%%%
H1 = L\(Z'*  (UQinv \ (LQinv \ (PQinv*(B/L')))));
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
      f = 10;
      g = 0;
      t = inf;
      tsdp = 0;
      return;
end
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%

      

%%%%%%%%%%%%%%
% Maximize lambda_max(H0 + t H1) over t using eigopt
% (Alternatively an SDP solver can be used for this purpose, but this seems much slower)
%%%%%%%%%%%%%%
parsin.H0 = H0;
parsin.H1 = H1;
parsin.tol = 10^-12;
parsin.gamma = -0.00000001;
parsin.minmax = 1;
intvalin.lb = -5*10^6;
intvalin.ub = 5*10^6;
      
tic;
[f,t,parsout] = eigopt('Hstructured_radii_inner_e',intvalin,parsin);
tsdp = toc;
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%


if (f >= 10)
      f = 10;
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
if modeQ
      
  if issparse(pars.J)
      W = pars.A - 1i*(omega+h)*speye(n);
  else
      W = pars.A - 1i*(omega+h)*eye(n);
  end
  Z = W\B;
      
else
      
  Qinv = pars.Qinv;
  if ~modeQs
      Z = Qinv*(((J-R)-(omega+h)*1i*Qinv)\B);
  else
      pars.s = 1i*(omega+h);
      pars.B = B;
      pars.Qinv = Qinv;
      pars.sq = 0;
      [Z,~] = feval(mylu,(J-R)-(omega+h)*1i*Qinv,pars);
      
  end
end


% FORM H0 and THE CHOLESKY FACTOR
H0 = B'* (UQinv \ (LQinv \ (PQinv*Z)));
H0 = H0'*H0;
H0 = (H0 + H0')/2;
% L = chol(H0,'lower');
 Rc = chol(H0);
 L = Rc';
      
H0 = L\(L'\eye(m));
H0 = (H0 + H0')/2;

% FORM H1
H1 = L\(Z'*  (UQinv \ (LQinv \ (PQinv*(B/L')))));
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
