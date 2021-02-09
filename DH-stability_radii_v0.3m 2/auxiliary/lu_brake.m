function [X,Y] = lu_brake(A,pars)



n = size(A,1);

s = pars.s;
B = pars.B;
Qinv = pars.Qinv;

m = size(B,2);


if ~isfield(pars,'sq')
    pars.sq = 1;
end


K = A(1:n/2,n/2+1:n);
M = A(1:n/2,1:n/2);


[L1,U1,P1] = lu(K);
[L2,U2,P2] = lu(-K-s*M);



%%%%%%%%%%%%%%%%
% (1) COMPUTE X1
%%%%%%%%%%%%%%%%
% keyboard
Bhat = [B(1:n/2,:); B(n/2+1:n,:)-s*B(1:n/2,:)];
X(n/2+1:n,:) = U2\(L2\(P2*Bhat(n/2+1:n,:)));
X(1:n/2,:) = U1\(L1\(P1*(Bhat(1:n/2,:) - M*X(n/2+1:n,:))));


%%%%%%%%%
% LAST but not the least
% Swap the variables
%%%%%%%%%

temp = X(1:n/2,:);
X(1:n/2,:) = X(n/2+1:n,:);
X(n/2+1:n,:) = temp;
X = Qinv*X;




%%%%%%%%%%%%%%%%
% (2) COMPUTE Y1
%%%%%%%%%%%%%%%%

if pars.sq
    Xhat = [X(1:n/2,:); X(n/2+1:n,:)-s*X(1:n/2,:)];
    Y(n/2+1:n,:) = U2\(L2\(P2*Xhat(n/2+1:n,:)));
    Y(1:n/2,:) = U1\(L1\(P1*(Xhat(1:n/2,:) - M*Y(n/2+1:n,:))));


    %%%%%%%%%
    % LAST but not the least
    % Swap the variables
    %%%%%%%%%

    temp = Y(1:n/2,:);
    Y(1:n/2,:) = Y(n/2+1:n,:);
    Y(n/2+1:n,:) = temp;
    Y = Qinv*Y;
else
    Y = zeros(n,m);
end


return;
