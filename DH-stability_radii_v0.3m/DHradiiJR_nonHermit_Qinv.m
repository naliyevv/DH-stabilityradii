function [f,z,info] = DHradiiJR_nonHermit_Qinv(J,Q,R,B,C,mode,initnum,initsub,intval,modeQ,modeQs,mylu)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% TASK:
% Computes the reciprocal of the stability radius for the
% dissipative-Hamiltonian (DH) system of the form
%                       x' =(J-R)Qx     (1)
% with respect to pertubations of one of J or R,
% using the structure-preserving subspace framework.
%
% The stability radii that are relevant are respectively defined by
%
%       inf { || Delta ||   |    (J + B (Delta) C - R)Q  has an eigenvalue on the imaginary axis  }
%       inf { || Delta ||   |    (J -  (R + B (Delta) C)Q  has an eigenvalue on the imaginary axis  }
%
% for given full-rank tall-skinny matrix B and short-fat matrix C.
%
% Both of the stability radii, when they are finite, has the characterization
%       1 / sup_{w in R} sigma_max( CQ (iwI - (J-R)Q)^{-1} B)
% where sigma_{max}(.) denotes the largest singular value of its matrix argument.
%
%
% ************************
% *** IMPORTANT REMARK ***
% ************************
% USE THIS ESPECIALLY IF Q^{-1} IS AVAILABLE AND SPARSE RATHER THAN Q.
% OTHERWISE IT MAY BE BETTER TO USE DHradiiJR_nonHermit.m
% ************************
% ************************
%
%
% CALL : [f,z,info] = DHradiiJR_nonHermit_Qinv(J,Q,R,B,C,mode,initnum,initsub,intval)
%
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
%     modeQ: to specify whether Q or Q^{-1} is available
%              modeQ = 0 - Q^{-1} is made available and passed through the input parameter Q.
%              modeQ = 1 - Q is made available.
%     (default - modeQ = 1)
%     modeQs: related to the formation of the right subspaces
%              (RS)       Col([ (i*omega*I - (J-R)Q)^{-1} B  (i*omega*I - (J-R)Q)^{-2} B ])
%             to interpolate at i*omega.
%               modeQs = 0 - Use Matlab's built-in linear system solvers
%               modeQs = 1 - User suplies the solver through the input parameter mylu
%     (default - modeQs = 0)
%     mylu: the name of the function written by the user to compute (RS), it should be
%           called as follows:
%               [X,Y] = mylu((i*omega*Qinv - (J-R)),pars)
%           the parameters, in particular B, should be passed through pars. At termination
%           X = Col([ (i*omega*I - (J-R)Q)^{-1} B   and   Y = Col([(i*omega*I - (J-R)Q)^{-2} B ])
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

if (nargin < 10)
    modeQ = 1;
end

if (nargin < 11)
    modeQs = 0;
end
                     


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form some of the system matrices and their factorizations that
% will be employed commonly
%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = speye(size(J));

if modeQ
    A = (J-R)*Q;
else
    Qinv = Q;
    [LQinv,UQinv,PQinv] = lu(Qinv);
end

% we assume D matrix is zero
D = zeros(size(C,1),size(B,2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FORM THE INITIAL SUBSPACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    intlength = intval.u - intval.l;
    h = intlength/(initnum-1);
    V = [];

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
            if modeQ
                maxsval(k) = max(svd(C*Q*((imag(eval(k))*E*1i  - (J-R)*Q)\B)));
            else
                maxsval(k) = max(svd(C*((imag(eval(k))*E*Qinv*1i  - (J-R))\B)));
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

        if modeQ
            [L,U] = lu(imag(eval(k))*1i*E-A);
            X1 = U\(L\B);
            Y1 = U\(L\X1);
        else
            if ~modeQs
                [L,U] = lu(imag(eval(k))*1i*E*Qinv-(J-R));
                X1 = Qinv*(U\(L\B));
                Y1 = Qinv*(U\(L\X1));
            else
                pars.s = 1i*imag(eval(k));
                pars.B = B;
                pars.Qinv = Qinv;
                [X1,Y1] = feval(mylu,imag(eval(k))*1i*E*Qinv-(J-R),pars);
            end
        end
        V = [V X1 Y1];
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%
% Form the bases for the right subspace and the left subspace
%%%%%%%%%%%%%%%%%%%%%
 V = full(V);
 
 % form orthonormal basis for the subspace V.
 [V,~] = qr(V,0);


 % construct the left subspace by means of Q and the right subspace V
 if modeQ
    W = Q * ( V/(V' * (Q * V)) );
 else
    VQV = V'*(UQinv\(LQinv\(PQinv*V)));
    W = UQinv\(LQinv \ (PQinv*(V/VQV)));
 end
 

 

%%%%%%%%%%%%%%%%%%%%%
% form the system matrices for the reduced PH system defined by the bases V and W
% corresponding to the right and left subspaces
%%%%%%%%%%%%%%%%%%%%%

 Jr = W' * J * W; 
 Rr = W' * R * W;
                 
 if modeQ
    Qr = V' * Q * V;
 else
    Qr = VQV;
 end
 
Br = W' * B;
Cr = C * W;
Ar = (Jr - Rr) * Qr;
Cr = Cr * Qr;
              
Br = full(Br);
Cr = full(Cr);
Ar = full(Ar);
                 

%%%%%%%%%%%%%
 % initializations for the main loop
 %%%%%%%%%%%%%
                            
 iter = 0;
              
 % prescribed tolerance
 tol = 10^-12;
              
 % maximum number of iterations
 maxit = 80;


 
%%%%%%%%%%%%%%%%%%%
% solve the reduced problem with the boyd-balakrishnan algorithm
% that is maximize sigma_max( Cr Qr (iwI - (Jr-Rr)Qr) Br) over w
% f is the globally maximal value, z is the global maximizer
%%%%%%%%%%%%%%%%%%%
 
 parssys = ss( Ar, Br, Cr, D );
 [f,z] = getPeakGain(parssys,0.05*tol);
          
                 
 zold = z/2;
 
 % fold will keep the computed value of the unstructured stability radius 
 % of the previous reduced system
 fold = f/2;
 
 
 
     
 %%%%%%%%%%%
 % main loop
 %%%%%%%%%%%
 while (abs(z-zold) > tol*abs(z) && abs(f-fold) > tol*f && iter < maxit)  

       %%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%
       % Expand the subspaces so that Hermite interpolation is achieved at i*z
       % where z is the global maximizer of sigma_max( Cr Qr (iwI - (Jr-Rr)Qr) Br) over w
       %%%%%%%%%%%%%%%%%%%%%%
        
        fold = f;
        zold = z;
                 
                 
        % solve new linear systems to expand the subspaces
        if modeQ
            [L,U] = lu(z*1i*E-A);
            Xz = U\(L\B);
            Yz = U\(L\Xz);
        else
                 
            if ~modeQs
                 [L,U] = lu(z*1i*E*Qinv-(J-R));
                 
                 Xz = Qinv*(U\(L\B));
                 Yz = Qinv*(U\(L\Xz));
            else
                 pars.s = 1i*z;
                 pars.B = B;
                 pars.Qinv = Qinv;
                 [Xz,Yz] = feval(mylu,z*1i*E*Qinv-(J-R),pars);
            end
                 
        end

        % expand the right subspace V
        V = [V Xz Yz];
        V = full(V);
   
        % form orthonormal basis for V
        [V,~] = qr(V,0);
        
        % construct the left subspace by means of Q and the right subspace V
        if modeQ
            W = Q * ( V/(V' * (Q * V)) ) ;
        else
            VQV = V'*(UQinv\(LQinv\(PQinv*V)));
            W = UQinv\(LQinv \ (PQinv*(V/VQV)));
        end
                         
      %%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%
 

        %%%%%%%%%%%%%%%%%%%%%
        % form the system matrices for the reduced PH system defined by the bases V and W
        % corresponding to the right and left subspaces
        %%%%%%%%%%%%%%%%%%%%%
        
        % reduced matrices
        Jr = W' * J * W;
                         
        if modeQ
            Qr = V' * Q * V;
        else
            Qr = VQV;
        end
        Rr = W' * R * W;
        Br = W' * B;
        Cr = C * W;

 
        Ar = (Jr - Rr) * Qr;
        Cr = Cr * Qr;
                     
        Br = full(Br);
        Cr = full(Cr);
        Ar = full(Ar);
                         
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%
                         
                         
        %%%%%%%%%%%%%%%%%%%
        % solve the reduced problem with the boyd-balakrishnan algorithm
        % that is maximize sigma_max( Cr Qr (iwI - (Jr-Rr)Qr) Br) over w
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
% its reciprocal is the unstrcutured distance to instability
%%%%%%
%%%%%%
                         
 
t2 = cputime;
ctime = t2-t1;
                         

info.iteration = iter;
info.time = ctime;
info.subspacedim = size(Ar,1);

return


 
