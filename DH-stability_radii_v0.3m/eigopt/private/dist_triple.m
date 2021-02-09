function [f, g] = dist_triple(gamma,pars)
% Emre Mengi (Modified August 19, 2011)
%
% call: [f, g] = dist_triple(gamma,pars)


A = pars.A;
lambda = pars.lambda;





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

	f = -S(3);
					
	g(1,1) = -real(um(1:n)'*vm(n+1:2*n));
	g(4,1) = +imag(um(1:n)'*vm(n+1:2*n));

	g(2,1) = -real(um(1:n)'*vm(2*n+1:3*n));
	g(5,1) = +imag(um(1:n)'*vm(2*n+1:3*n));

	g(3,1) = -real(um(n+1:2*n)'*vm(2*n+1:3*n));
	g(6,1) = imag(um(n+1:2*n)'*vm(2*n+1:3*n));
else
	[eigvec,d] = eigs(B'*B,3,'sm');
				  
	vm = eigvec(1:3*n,3);
	vm = vm/norm(vm);
				  
	f = -sqrt(real(d(3,3)));
				  
	um = B*vm;
	um = um/norm(um);
				  
	g(1,1) = -real(um(1:n)'*vm(n+1:2*n));
	g(4,1) = +imag(um(1:n)'*vm(n+1:2*n));
			   
	g(2,1) = -real(um(1:n)'*vm(2*n+1:3*n));
	g(5,1) = +imag(um(1:n)'*vm(2*n+1:3*n));
			   
	g(3,1) = -real(um(n+1:2*n)'*vm(2*n+1:3*n));
	g(6,1) = imag(um(n+1:2*n)'*vm(2*n+1:3*n));
end



return;