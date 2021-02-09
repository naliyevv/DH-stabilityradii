function [f,g,t,tsdp,tchol] = Hstructured_radii_inner_e(t,pars)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% Auxiliary routine called by Hstructured_radii_e, Hstructured_radii_e_invQs
%
% TASK:
% Computes lambda_min(pars.H0 + t pars.H1) and its derivative at a given t,
% returns the computed values inside f and g.


[eigvecs,eigvals] = eig(pars.H0 + t*pars.H1);

[f,indx] = min(real(diag(eigvals)));

g = real(eigvecs(:,indx)'*pars.H1*eigvecs(:,indx));

t = 0;
tsdp = 0;
tchol = 0;

return
