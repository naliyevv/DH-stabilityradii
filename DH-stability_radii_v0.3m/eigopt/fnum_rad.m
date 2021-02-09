function [f, g] = fnum_rad(theta,pars)
% Emre Mengi (Modified on August 19, 2011)
%
% call: function [f, g] = fnum_rad(theta,pars)
% task:
%		calculates f = lambda_max((A*e^(i*theta) + A'*e^(-i*theta))/2) 
%		and its derivative g
% note:
%		the input matrix A must be passed through pars.A



A = pars.A;
[n,m] = size(A);


if (n ~= m)
	error('the input matrix must be square');
end

H = (A*exp(i*theta) + A'*exp(-i*theta))/2;


if (n <= 250)
	 [V,D] = eig(H);

	 D1 = diag(real(D));	 
	 [f,ind] = max(D1);

	 g = real(V(:,ind)'*((i*A*exp(i*theta) - i*A'*exp(-i*theta))/2)*V(:,ind));
else
	 [eigvec,d] = eigs(H,1,'lr');	 
						   
	 vm = eigvec(:,1);
	 vm = vm/norm(vm);
						   
	 f = d(1,1);
	 g = real(vm'*((i*A*exp(i*theta) - i*A'*exp(-i*theta))/2)*vm);
end
	
					
					

return;