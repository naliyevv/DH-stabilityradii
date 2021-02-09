function [xvals,fvals] = plot_quad(xk,fk,gk,gamma)
% Emre Mengi (June 16, 2013)
%

xvals = [];
fvals = [];

for x = 0:0.003:2*pi
	f = fk + gk'*(x-xk) + (gamma/2)*(norm(x-xk))^2;

	xvals = [xvals x];
	fvals = [fvals f];
end


return;