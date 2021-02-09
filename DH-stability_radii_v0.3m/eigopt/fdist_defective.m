function [f, g] = fdist_defective(lambda,pars)
% Emre Mengi (Modified on August 19, 2011)
%
% call: function [f, g] = fdist_defective(lambda,pars)
% task:
%		calculates 
%				f = sup_gamma 
%						sigma_(2n-1)([A - (lambda(1)+i*lambda(2))*I		gamma*I;
%													0		A - (lambda(1)+i*lambda(2))*I])
%		and its gradient g
% note:
%		the input matrix A must be passed through pars.A





[n,n2] = size(pars.A);

if (n ~= n2)
	error('A must be square');
end

%%%%%%%%%%%%%
% (1) First find optimal gamma
%%%%%%%%%%%%%
pars.lambda = lambda(1) + i*lambda(2);
pars.nvar = 2;
pars.fgname = 'dist_defective';


options.prtlevel = 0;
options.x0 = [1 -1]';
options.normtol = 10^-10;

[gamma,f] = bfgs(pars,options);






% (2) Now evaluate the supremum function and
%	  its derivative at the optimal gamma
A = pars.A;

[n,n2] = size(A);

B = [A-pars.lambda*eye(n)   (gamma(1)+i*gamma(2))*eye(n);
zeros(n)		A-pars.lambda*eye(n)];


if (n <= 150)
	[U,S,V] = svd(B);

	S = diag(S);
	[S,ind] = sort(S);

	vm = V(:,ind(2));
	vm = vm/norm(vm);
	um = U(:,ind(2));
	um = um/norm(um);

	f = S(2);
	g(1,1) = -real(um'*vm);
	g(2,1) = imag(um'*vm);
else
	[eigvec,d] = eigs(B'*B,2,'sm');

	vm = eigvec(1:2*n,2);
	vm = vm/norm(vm);

	f = sqrt(real(d(2,2)));

	um = B*vm;
	um = um/norm(um);

	g(1,1) = -real(um'*vm);
	g(2,1) = imag(um'*vm);
end





return;