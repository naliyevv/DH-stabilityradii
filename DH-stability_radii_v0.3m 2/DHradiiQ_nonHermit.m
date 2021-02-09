function [f,z,info] = DHradiiQ_nonHermit(J,Q,R,B,C,mode,initnum,initsub,intval)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% TASK:
% Computes the reciprocal of the stability radius for the
% dissipative-Hamiltonian (DH) system of the form
%                       x' =(J-R)Qx     (1)
% with respect to pertubations of Q,
% using the structure-preserving subspace framework.
%
% The relevant stability radius is formally defined by
%
%       inf { || Delta ||   |    (J - R) (Q + B Delta C)  has an eigenvalue on the imaginary axis  }
%
% for given full-rank tall-skinny matrix B and short-fat matrix C.
%
% The stability radius, when it is finite, has the characterization
%       1 / sup_{w in R} sigma_max( C (iwI - (J-R)Q)^{-1} (J-R)B)
% where sigma_{max}(.) denotes the largest singular value of its matrix argument.
%
%
%
% CALL : [f,z,info] = DHradiiQ_nonHermit(J,Q,R,B,C,mode,initnum,initsub,intval)
%
% INPUT:
%     J: skew-symmetric matrix for the energy flux of the system
%     Q: symmetric, positive definite matrix for the total energy of the system
%     R: symmetric positive semi-definite dissipation matrix of the system
%     B: restriction matrix from the left-hand side (must be full-rank)
%     C: restriction matrix from the right-hand side (must be full-rank)
%     mode: this is for the selction of the initial interpolation points.
%           mode = 0 : randomly chosen initnum of points iw where w in [intval.l, intval.u]
%           mode = 1 : [intval.l, intval.u] is first split into initnum of equally-spaced points,
%                      w_1, ..., w_initnum, then for each one of i*w_1, ..., i*w_initnum, the
%                      nearest eigenvalue of (J-R)Q is computed. The imaginary parts of these eigenvalues
%                      form the maximal set of initial interpolation points.
%           mode = 2 : [intval.l, intval.u] is split into initnum of equally-spaced points w_1, ..., w_initnum,
%                      the maximal set of initial interpolation points are given by i*w_1, ..., i*w_initnum
%     (default - mode = 2)
%     initnum: maximal number of initial interpolation points
%     (default - initnum = 10)
%     initsub: among initnum of initial interpolation points initsub of them are selected
%              based on which of these have larger sigma_max( C (iwI - (J-R)Q)^{-1} (J-R)B)
%     (default - initsub = round(initnum/4))
%     intval: specifies the interval where initial interpolation is to be performed
%             with the fields
%               intval.l, intval.u
%               the interval for initial interpolation is [i*intval.l, i*intval.u]
%     (default - [intval.l, intval.u] = [-2, 0 ])
%
%
% OUTPUT:
%       f: computed reciprocal of the unstructured stability radius of the DH system, i.e.,
%           f normally should be equal to sup_{w in R} sigma_max( C (iwI - (J-R)Q)^{-1} (J-R)B)
%           up to tolerances.
%       z: normally should be equal to argmax_{w in R} sigma_max( C (iwI - (J-R)Q)^{-1} (J-R)B)
%    info: structure containing information regarding the progress of the algorithm
%           info.iteration - total number of iterations performed to reach prescribed accuracy
%           info.time: total time consumed by the algorithm to reach prescribed accuracy
%           info.subspacedim: the dimension of the subspaces for projection at termination


warning off;
t1 = cputime;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form the port-Hamiltonian system associated with the transfer function H(s) = C (sI - (J-R)Q)^{-1} (J-R)B
% whose H-infinity norm is to be computed
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if issparse(J)
    E = speye(size(J));
else
    E = eye(size(J));
end

% D matrix is always zero
D = zeros(size(C,1),size(B,2));

A = (J-R) * Q;
JRB = (J-R) * B;
Ct = C';
 


%%%%%%%%%%%%%%%%%%%%%%
% Set the default values of the parameters
%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 6)
    mode = 2;
end

if (nargin < 7)
   initnum = 10;
end

if (nargin < 8)
    initsub = round(initnum/4);
end

if (nargin < 9)
    intval.l = -2;
    intval.u = 0;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FORM THE INITIAL SUBSPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    intlength = intval.u - intval.l;
    h = intlength/(initnum-1);
    W = [];

    for k = 1:initnum
        if mode
            w = intval.l + (k-1)*h;
            if mode == 1
                eval(k) = eigs(A,1,w*1i);
            else
                eval(k) = w*1i;
            end
        else
            eval(k) = (intval.l + rand*intlength)*1i;
        end

        if (initnum > initsub)
            if issparse(J)
                MM = C*((imag(eval(k))*E*1i  - A)\JRB);
                maxsval(k) = max(svd(full(MM)));
            else
                maxsval(k) = max(svd(C*((imag(eval(k))*E*1i  - A)\JRB)));
            end
        end
    end

    if (initnum > initsub)
        [maxsval,indx] = sort(maxsval,'descend');
        eval = eval(indx);
    end

    %%%%%%%%%%%%%%%%%%%%%%
    % at this point, the initial interpolation points are determined, namely
    % i*eval(1), i*eval(2), ..., i*eval(initsub)
    %%%%%%%%%%%%%%%%%%%%%%


    for k = 1:initsub
        info.eval(k) = eval(k);

        if (initnum > initsub)
            info.maxsval(k) = maxsval(k);
        else
            info.maxsval(k) = -1;
        end


        [L,U] = lu(imag(eval(k))*1i*E-A);
        Lt = L';
        Ut = U';
        X1 = Lt\(Ut\Ct);
        Y1 = Lt\(Ut\X1);

        W = [W X1 Y1];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
% Form the bases for the left subspace and the right subspace
%%%%%%%%%%%%%%%%%%%%%
 
 W = full(W);
 
 % form orthonormal basis for the subspace W.
 [W,~] = qr(W,0);
 
 % construct the right subspace by means of (J-R) and the left subspace W
 V = (J-R)' * ( W/(W' * ((J-R)' * W)) ) ;
 

                         
 


 %%%%%%%%%%%%%%%%%%%%%
 % form the system matrices for the reduced PH system defined by the bases V and W
 % corresponding to the right and left subspaces
 %%%%%%%%%%%%%%%%%%%%%
 Jr = W' * J * W; 
 Rr = W' * R * W;
 Qr = V' * Q * V;
 
              
Br = V' * B;
Cr = C * V;
Ar = (Jr - Rr) * Qr;
Br = (Jr - Rr) * Br;
                         
Br = full(Br);
Cr = full(Cr);
Ar = full(Ar);


 %%%%%%%%%%%%%
 % initializations for the main loop
 %%%%%%%%%%%%%

 % fold will keep the computed value of the unstructured stability radius of the
 % previous reduced system


 fold = 0;
                            
 iter = 0;
              
 % prescribed tolerance
 tol = 10^-12;
              
 % maximum number of iterations
 maxit = 80;
                         
                         
%%%%%%%%%%%%%%%%%%%
% solve the reduced problem with the boyd-balakrishnan algorithm
% that is maximize sigma_max( Cr (iwI - (Jr-Rr)Qr) (Jr - Rr)Br) over w
% f is the globally maximal value, z is the global maximizer
%%%%%%%%%%%%%%%%%%%
 
 parssys = ss( Ar, Br, Cr, D );
 [f,z] = getPeakGain(parssys,0.05*tol);

 zold = z-1;
 
 
 
 %%%%%%%%%%%
 % main loop
 %%%%%%%%%%%
 while (abs(z-zold) > tol*abs(z) && abs(f-fold) > tol*f && iter < maxit)


        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        % Expand the subspaces so that Hermite interpolation is achieved at i*z
        % where z is the global maximizer of sigma_max( Cr  (iwI - (Jr-Rr)Qr) (Jr-Rr)Br ) over w
        %%%%%%%%%%%%%%%%%%%%%%
                         
        fold = f;
        zold = z;
        
                         
        % solve new linear systems to expand the subspaces
                     
        [L,U] = lu(z*1i*E-A);
        Lt = L';
        Ut = U';
        Xz = Lt\(Ut\Ct);
        Yz = Lt\(Ut\Xz);

        % expand the left subspace
        w = [Xz Yz];
        W = [W w];
        W = full(W);
                 
        % form orthonormal basis for W
        [W,~] = qr(W,0);
                 
        % construct the right subspace by means of Q and the left subspace W
        V = (J-R)' * ( W/(W' * ((J-R)' * W)) ) ;
                                
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
 

 

       %%%%%%%%%%%%%%%%%%%%%
       % form the system matrices for the reduced PH system defined by the bases V and W
       % corresponding to the right and left subspaces
       %%%%%%%%%%%%%%%%%%%%%
                                
        
        % reduced matrices
        Jr = W' * J * W; 
        Qr = V' * Q * V;
        Rr = W' * R * W;
        Br = V' * B;
        Cr = C * V;

 
        Ar = (Jr - Rr) * Qr;
        Br = (Jr - Rr) * Br;
                     
        Br = full(Br);
        Cr = full(Cr);
        Ar = full(Ar);

        %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
                     
        %%%%%%%%%%%%%%%%%%%
        % solve the reduced problem with the boyd-balakrishnan algorithm
        % that is maximize sigma_max( Cr (iwI - (Jr-Rr)Qr) (Jr-Rr)Br) over w
        % f is the globally maximal value, z is the global maximizer
        %%%%%%%%%%%%%%%%%%%
        parssys = ss( Ar, Br, Cr, D );
        [f,z] = getPeakGain(parssys,0.05*tol);
        
        iter = iter + 1;
          
 end

                                
%%%%%%
%%%%%%
% NOTE:
% the computed f is the H-infinity norm of the PH system
% its reciprocal is the non-Hermitian stability radius
%%%%%%
%%%%%%
 
                                
t2 = cputime;
ctime = t2-t1;
                                

info.iteration = iter;
info.time = ctime;
info.subspacedim = size(Ar,1);

return


 
