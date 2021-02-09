function [x, f, g, H] = bfgs1run(x0, pars, options)
% Copyright Michael Overton
% make a single run of BFGS from one starting point
% intended to be called by bfgs only
% reference: Nocedal and Wright
fgname = pars.fgname;
normtol = options.normtol;
fvalquit = options.fvalquit;
cpufinish = cputime + options.cpumax;
maxit = options.maxit;
prtlevel = options.prtlevel;
strongwolfe = options.strongwolfe;
wolfe1 = options.wolfe1;
wolfe2 = options.wolfe2;
x = x0;
% keyboard
[f,g] = feval(fgname, x, pars);
gnorm = norm(g);
H0 = options.H0;
if isempty(H0)
    H = eye(length(x0));
    scaleinit = 1; % so H is scaled before the first update
else
    H = H0;
    scaleinit = 0; % no scaling if H0 is provided
end
if f == inf % better not to generate an error return
    if prtlevel > 0
        fprintf('bfgs: f is infinite at initial iterate\n')
    end
    return
elseif isnan(f)
    if prtlevel > 0
        fprintf('bfgs: f is nan at initial iterate\n')
    end
    return
elseif gnorm < normtol
    if prtlevel > 0
        fprintf('bfgs: tolerance on gradient satisfied at initial iterate\n')
    end
    return
elseif f < fvalquit
    if prtlevel > 0
        fprintf('bfgs: below target objective at initial iterate\n')
    end
    return
end 
for iter = 1:maxit
%	keyboard
    p = -H*g; % H approximates the inverse Hessian
    gtp = g'*p;
    if gtp >= 0
       if prtlevel > 0
          fprintf('bfgs: not a descent direction, quit at iteration %d, f = %e, gnorm = %5.1e\n',...
              iter, f, gnorm)
       end
       return
    end
    if strongwolfe % strong Wolfe line search is optional
        [alpha, x, f, gnew, fail] = ...
            linesch_sw(x, f, gtp, p, pars, wolfe1, wolfe2, prtlevel);
    else % weak Wolfe line search is the default
		fold = f;
        [alpha, x, f, gnew, fail] = ...
            linesch_ww(x, f, gtp, p, pars, wolfe1, wolfe2, prtlevel);
    end
    if prtlevel > 1
        fprintf('bfgs: iter %d: step = %e, f = %e, gnorm = %5.1e\n', iter, alpha, f, gnorm)
    end
    if fail == 1 % Wolfe conditions not both satisfied, quit
        if prtlevel > 0
           fprintf('bfgs: quit at iteration %d, f = %e, gnorm = %5.1e\n',...
               iter, f, gnorm)
        end
        return
    elseif fail == -1 % function apparently unbounded below
        if prtlevel > 0
           fprintf('bfgs: f may be unbounded below, quit at iteration %d, f = %e\n', iter, f)
        end
        return
    end
    gnorm = norm(gnew);
    if gnorm <= normtol 
        if prtlevel > 0
            fprintf('bfgs: gradient norm below tolerance, quit at iteration %d, f = %e\n', iter, f')
        end
        g = gnew;
        return
    end
	if abs(f - fold) < normtol
		if prtlevel > 0
			fprintf('bfgs: function values too close %d, f = %e %e\n', iter, f,fold)
        end
		g = gnew;
		return
	end
    if f < fvalquit
        if prtlevel > 0
            fprintf('bfgs: reached target objective, quit at iteration %d \n', iter)
        end
        return
    end
    if cputime > cpufinish
        if prtlevel > 0
            fprintf('bfgs: quit since CPU time limit exceeded\n')
        end
        return
    end
    s = alpha*p;
    y = gnew - g;
    sty = s'*y;    % successful line search ensures this is positive
    if sty > 0     % perform rank two BFGS update
        if iter == 1 & scaleinit % Nocedal and Wright recommend scaling I before first update
            H = (sty/(y'*y))*H; % equivalently, replace H on the right by I
        end
        rho = 1/sty;
        rhoHyst = rho*(H*y)*s';                                       % M = I - rho*s*y';
        H = H - rhoHyst' - rhoHyst + rho*s*(y'*rhoHyst) + rho*s*s';   % H = M*H*M' + rho*s*s';
        H = 0.5*(H + H');   % may not be symmetric because of rounding
    else % should not happen unless line search fails
        if prtlevel > 0
            fprintf('bfgs: *** sty <= 0, skipping BFGS update at iteration %d \n', iter) 
        end
    end
    g = gnew;
end % for loop
if prtlevel > 0
    fprintf('bfgs: %d iterations reached, f = %e, gnorm = %5.1e\n', maxit, f, gnorm)
end