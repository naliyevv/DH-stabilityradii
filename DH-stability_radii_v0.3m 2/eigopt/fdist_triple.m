function [f, g] = fdist_triple(lambda,pars)
% Emre Mengi (Modified on August 19, 2011)
%
% call: function [f, g] = fdist_triple(lambda,pars)
% task:
%		calculates 
%				f = sup_gamma 
%						sigma_(3n-2)([A - (lambda(1)+i*lambda(2))*I		gamma_12*I			gamma_13*I;
%													0		A - (lambda(1)+i*lambda(2))*I	gamma_23*I;
%													0			0							A - (lambda(1)+i*lambda(2)*I])
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
pars.nvar = 6;
pars.fgname = 'dist_triple';


options.prtlevel = 0;
options.x0 = [1 1 1 -1 -1 -1]';
options.normtol = 10^-10;

[gamma,f] = bfgs(pars,options);






% (2) Now evaluate the supremum function and
%	  its derivative at the optimal gamma
A = pars.A;

[n,n2] = size(A);
Zn = zeros(n);

B = [A-pars.lambda*eye(n)   (gamma(1)+i*gamma(4))*eye(n)		(gamma(2)+i*gamma(5))*eye(n);
				Zn				A-pars.lambda*eye(n)			(gamma(3)+i*gamma(6))*eye(n);
				Zn				Zn								A-pars.lambda*eye(n)];


if (n <= 150)
	[U,S,V] = svd(B);

	S = diag(S);
	[S,ind] = sort(S);

	vm = V(:,ind(3));
	vm = vm/norm(vm);
	um = U(:,ind(3));
	um = um/norm(um);

	f = S(3);
	g(1,1) = -real(um'*vm);
	g(2,1) = imag(um'*vm);
else
	[eigvec,d] = eigs(B'*B,3,'sm');

	vm = eigvec(1:3*n,3);
	vm = vm/norm(vm);

	f = sqrt(real(d(3,3)));

	um = B*vm;
	um = um/norm(um);

	g(1,1) = -real(um'*vm);
	g(2,1) = imag(um'*vm);
end





return;