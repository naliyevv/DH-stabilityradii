function [f, g] = dist_defective(gamma,pars)
% Emre Mengi (Modified August 19, 2011)
%
% call: [f, g] = dist_defective(gamma,pars)


A = pars.A;
lambda = pars.lambda;





[n,n2] = size(A);

B = [A-lambda*eye(n)   (gamma(1)+i*gamma(2))*eye(n);
	 zeros(n)		A-lambda*eye(n)];





if (n <= 150)
	[U,S,V] = svd(B);
				  
	S = diag(S);
	[S,ind] = sort(S);

	vm = V(:,ind(2));
	vm = vm/norm(vm);
	um = U(:,ind(2));
	um = um/norm(um);

	f = -S(2);
					
	g(1,1) = -real(um(1:n)'*vm(n+1:2*n));
	g(2,1) = +imag(um(1:n)'*vm(n+1:2*n));
else
	[eigvec,d] = eigs(B'*B,2,'sm');
				  
	vm = eigvec(1:2*n,2);
	vm = vm/norm(vm);
				  
	f = -sqrt(real(d(2,2)));
				  
	um = B*vm;
	um = um/norm(um);
				  
	g(1,1) = -real(um(1:n)'*vm(n+1:2*n));
	g(2,1) = +imag(um(1:n)'*vm(n+1:2*n));
end



return;