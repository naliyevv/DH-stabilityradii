function [f,g] = H_infinity(s,pars)
% Emre Mengi (Modified on August 19, 2011)
%
% call: function [f, g] = H_infinity(s,pars)
% task:
%		calculates f = sigma_max(C*(s*I - A)^-1*B) and its
%		derivative g
% note:
%		the input matrices A,B,C must be passed through 
%		pars.A, pars.B, pars.C.


A = pars.A;
B = pars.B;
C = pars.C;
D = pars.D;

[m,n]=size(A);

if (n ~= m)
	error('A must be square');
end

[nr,nc] = size(B);

if (nr ~= n)
	error('A and B must have same number of rows');
end

[nr,nc] = size(C);

if (nc ~= n)
	error('A and C must have same number of columns');
end





H = C*((s*i*eye(m)-A)^-1)*B + D;

[U,S,V] = svd(H);

f = S(1,1);

g = imag(U(:,1)'*C*((s*i*eye(m)-A)^-2)*B*V(:,1));



return



