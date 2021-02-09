function [f,z,info] = DHradiiJR_Hermit(J,Q,R,B,mode,initnum,initsub,intval)
%
% Copyright: N. Aliyev, V. Mehrmann, E. Mengi
%
% TASK:
% Computes the structured stability radius for the
% dissipative-Hamiltonian (DH) system of the form
%                       x' =(J-R)Qx     (1)
% with respect to Hermitian pertubations of R,
% using the structure-preserving subspace framework.
%
% The structured stability radius is defined by
%
%    inf { || Delta || | Delta = Delta' and (J -  (R + B Delta B'))Q has an eigenvalue on the imaginary axis  }
%
% for a given full-rank tall-skinny matrix B.
%
% The subspace framework is based on the eigenvalue optimization characterization in [1, Theorem 3.2].
%
%
% CALL : [f,z,info] = DHradiiJR_Hermit(J,Q,R,B,mode,initnum,initsub,intval)
%
%
% INPUT:
%     J: skew-symmetric matrix for the energy flux of the system
%     Q: symmetric, positive definite matrix for the total energy of the system
%     R: symmetric positive semi-definite dissipation matrix of the system
%     B: the restriction matrix (must be full-rank)
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
%              based on which of these make the objective eigenvalue (to be minimized) smaller.
%     (default - initsub = round(initnum/4))
%     intval: specifies the interval where initial interpolation is to be performed
%             with the fields
%               intval.l, intval.u
%               the interval for initial interpolation is [i*intval.l, i*intval.u]
%     (default - [intval.l, intval.u] = [-2, 0 ])
%
%
% OUTPUT:
%       f: computed structured stability radius of the DH system, i.e.,
%       z: the real number such that i*z becomes an eigenvalue under the smallest Hermitian perturbation
%          possible of norm equal to the structured stability radius.
%    info: structure containing information regarding the progress of the algorithm
%           info.iteration - total number of iterations performed to reach prescribed accuracy
%           info.time: total time consumed by the algorithm to reach prescribed accuracy
%           info.subspacedim: the dimension of the subspaces for projection at termination
%
%
% REFERENCES:
% [1] - Computation of Stability Radii for Large-Scale Port-Hamiltonian Systems. Tech. Report.
%       Nicat Aliyev, Volker Mehrmann and Emre Mengi.

warning off;
t1 = cputime;


%%%%%%%%%%%%%%%%%
% Some Initializations
%%%%%%%%%%%%%%%%%

pars.J = J;
pars.Q = Q;
pars.R = R;
pars.B = B;
A = (J-R)*Q;
pars.A = A;
pars.tol = 10^-7;
pars.gamma = -20;
intval.lb = intval.l;
intval.ub = intval.u;

if issparse(J)
    In = speye(size(J));
else
    In = eye(size(J));
end


%%%%%%%%%%%%%%%%%%%%%%
% Set the default values of the parameters
%%%%%%%%%%%%%%%%%%%%%%
if (nargin < 5)
    mode = 2;
end

if (nargin < 6)
   initnum = 10;
end

if (nargin < 7)
   initsub = round(initnum/4);
end

if (nargin < 8)
   intval.l = -2;
   intval.u = 0;
end




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
            mineigval(k) = Hstructured_radii_e(imag(eval(k)),pars);
        end
    end

    if (initnum > initsub)
        [mineigval,indx] = sort(mineigval,'ascend');
        eval = eval(indx);
    end

    %%%%%%%%%%%%%%%%%%%%%%
    % at this point, the initial interpolation points are determined, namely
    % i*eval(1), i*eval(2), ..., i*eval(initsub)
    %%%%%%%%%%%%%%%%%%%%%%


    for k = 1:initsub
        info.eval(k) = eval(k);


        if (initnum > initsub)
            info.mineigval(k) = mineigval(k);
        else
            info.mineigval(k) = -1;
        end
    
        [L,U] = lu(A - imag(eval(k))*1i*In);
        X1 = U\(L\B);
        Y1 = U\(L\X1);
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
 W = Q * ( V/(V' * (Q * V)) ) ;
 

 
%%%%%%%%%%%%%%%%%%%%%
% form the matrices for the reduced problem defined by the bases V and W
% corresponding to the right and left subspaces
%%%%%%%%%%%%%%%%%%%%%
 
 Jr = W' * J * W; 
 Rr = W' * R * W;
 Qr = V' * Q * V;
 Br = W' * B;

              



              
            
%%%%%%%%%%%%%%%%%%%
% solve the reduced problem using eigopt
% f is the globally minimal value of the eigenvalue objective, z is the global minimizer
%%%%%%%%%%%%%%%%%%%

pars.J = full(Jr);
pars.Q = full(Qr);
pars.R = full(Rr);
pars.A = (pars.J - pars.R)*pars.Q;
pars.B = full(Br);
              
[f,z,parsout] = eigopt('Hstructured_radii_e',intval,pars);
              
info.tsdp = parsout.tsdp;
info.tchol = parsout.tchol;
              
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
              


              
%%%%%%%%%%%%%
% initializations for the main loop
%%%%%%%%%%%%%

  % fold will keep the computed value of the unstructured stability radius of the
  % previous reduced system
  fold = 0;
              
  zold = z - 1;
 
 % initial iteration
 iter = 0;
 
 % prescribed tolerance
 tol = 10^-6;
 
 % maximum number of iterations
 maxit = 100;
 
 
              
%%%%%%%%%%%
% main loop
%%%%%%%%%%%
 while (abs(z-zold)>tol*abs(z) && abs(f-fold)>tol*f && iter < maxit)
              

       %%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%
       % Expand the subspaces so that Hermite interpolation is achieved at i*z
       % where z is the global maximizer of sigma_max( Cr Qr (iwI - (Jr-Rr)Qr) Br) over w
       %%%%%%%%%%%%%%%%%%%%%%
        
        fold = f;
        zold = z;
        
              
        % solve linear systems to expand the subspaces
        [L,U] = lu(A - z*1i*In);
        Xz = U\(L\B);
        Yz = U\(L\Xz);
              
        % expand the right subspace V
        V = [V Xz Yz];
        V = full(V);
   
        % form orthonormal basis for V
        [V,~] = qr(V,0);
        
        % construct the left subspace by means of Q and the right subspace V
        W = Q * ( V/(V' * (Q * V)) );
 
       %%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%

        
       %%%%%%%%%%%%%%%%%%%%%
       % form the matrices for the reduced problem defined by the bases V and W
       % corresponding to the right and left subspaces
       %%%%%%%%%%%%%%%%%%%%%
                     
        % reduced matrices
        Jr = W' * J * W; 
        Qr = V' * Q * V;
        Rr = W' * R * W;
        Br = W' * B;

        pars.J = full(Jr);
        pars.Q = full(Qr);
        pars.R = full(Rr);
        pars.A = (pars.J - pars.R)*pars.Q;
        pars.B = full(Br);
                     
                     

                     
                     
        %%%%%%%%%%%%%%%%%%%
        % solve the reduced problem using eigopt
        % f is the globally minimal value of the eigenvalue objective, z is the global minimizer
        %%%%%%%%%%%%%%%%%%%
        [f,z,parsout] = eigopt('Hstructured_radii_e',intval,pars);
        
        info.tsdp = info.tsdp + parsout.tsdp;
        info.tchol = info.tchol + parsout.tchol;
                     
        iter = iter + 1;
     
 end
                     

f = sqrt(f);
 
t2 = cputime;
ctime = t2-t1;
        
                     
info.iteration = iter;
info.time = ctime;
info.subspacedim = size(Jr,1);

return


 
