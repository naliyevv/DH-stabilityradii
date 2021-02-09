function [f, g] = fdist_uncont(lambda,pars)
% Emre Mengi (Modified on August 19, 2011)
%
% call: function [f, g] = fdist_uncont(lambda,pars)
% task:
%		calculates f = sigma_min([A - (lambda(1)+i*lambda(2))*I		B]) 
%		and its gradient g
% note:
%		the input matrices A,B must be passed through pars.A, pars.B


A = pars.A;
B = pars.B;

[n,n2] = size(A);
[nb,m] = size(B);

if (n ~= n2)
	error('A must be square');
end

if (n ~= nb)
	error('A and B must have same number of rows');
end


RM = [A-(lambda(1)+i*lambda(2))*eye(n)   B];

 if (n >= 150)
	[eigvec,d] = eigs(RM*RM',1,'sm');

	um = eigvec(1:n,1);
	um = um/norm(um);
				   
	f = sqrt(real(d(1,1)));
				  
	vm = RM'*um;
	vm = vm/norm(vm);

				  
	g(1,1) = -real(um'*vm(1:n));
	g(2,1) = imag(um'*vm(1:n));
 else
	[U,S,V] = svd(RM);
	
	f = real(S(n,n));
	g(1,1) = -real(U(:,n)'*V(1:n,n));
	g(2,1) = imag(U(:,n)'*V(1:n,n));
 end
				  
			 

return;