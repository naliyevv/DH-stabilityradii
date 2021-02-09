function [f, g] = fdist_instab(omega,pars)
% Emre Mengi (Modified on August 19, 2011)
%
% call: function [f, g] = fdist_instab(omega,pars)
% task:
%		calculates f = sigma_min(A - omega*i*I) and its
%		derivative g
% note:
%		the input matrix A must be passed through pars.A



A = pars.A;
[n,m] = size(A);


if (n ~= m)
	error('the input matrix must be square');
end



if (n <= 150)
	[U,S,V] = svd(A - eye(n)*i*omega);

	f = S(n,n);
	g = imag(U(:,n)'*V(:,n));
else
	[eigvec,d] = eigs((A-omega*i*eye(n))'*(A - omega*i*eye(n)),1,'sm');

	vm = eigvec(1:n,1);
	vm = vm/norm(vm);

	f = sqrt(real(d(1,1)));

	um = (A-omega*i*eye(n))*vm;				  
	um = um/norm(um);				  

	g = imag(um'*vm);
end




return;