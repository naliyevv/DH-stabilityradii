function quadval = evalq(x,xk,fk,gk,gamma)
% Mustafa Kilic, Emre Mengi and E. Alper Yildirim
% (Modified version July 30, 2012)
%

quadval = fk + gk'*(x-xk) + (gamma/2)*(norm(x-xk))^2;

return;