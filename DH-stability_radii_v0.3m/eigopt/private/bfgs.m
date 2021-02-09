function [x, f, g, H] = bfgs(pars, options)
%Copyright Michael Overton
%BFGS The BFGS quasi-Newton minimization algorithm.
%   Calls:  [x, f, g, H] = bfgs(pars) 
%    and:   [x, f, g, H] = bfgs(pars, options)
% 
%   Input parameters
%    pars is a required struct, with two required fields
%      pars.nvar: the number of variables
%      pars.fgname: string giving the name of function (in single quotes) 
%         that returns the function and its gradient at a given input x, 
%         with call   [f,g] = fgtest(x,pars)  if pars.fgname is 'fgtest'.
%         Any data required to compute the function and gradient may be
%         encoded in other fields of pars.
%    options is an optional struct, with no required fields
%      options.x0: each column is a starting vector of variables
%          (default: generated randomly)
%      options.nrand: number of starting vectors to be generated randomly
%          in addition to those specified in options.x0, if any
%          (default: 10 if options.x0 is not specified, 0 if it is)
%      options.maxit: max number of iterations (default 100)
%      options.normtol: absolute tolerance on gradient norm (default: 1e-6)
%      options.fvalquit: quit if f drops below this value (default: -inf)
%      options.cpumax: quit if cpu time in secs exceeds this (default: inf)
%      options.H0: initial inverse Hessian approximation (default: identity)
%      options.strongwolfe: 0 for weak Wolfe line search (default)
%                           1 for strong Wolfe line search
%      options.wolfe1: first Wolfe line search parameter 
%          (ensuring sufficient decrease in function value, default: 1e-4)
%      options.wolfe2: second Wolfe line search parameter 
%          (ensuring algebraic increase (weak) or absolute decrease (strong)
%           in projected gradient, default: 0.9)
%          ("strong" Wolfe line search is not usually recommended for use with
%           BFGS; it is very complicated and bad if f is nonsmooth, but
%           it could have advantages in some cases)
%      options.prtlevel: one of 0 (no printing), 1 (minimal), 2 (verbose)
%          (default: 1)
%
%   Output parameters 
%    x: each column is an approximate minimizer, one for each starting vector
%    f: each entry is the function value for the corresponding column of x
%    g: each column is the gradient at the corresponding column of x
%    H: cell array of the final BFGS inverse Hessian approximation matrices
%     (except that H is a matrix if there is only one starting vector)
%
%   BFGS is normally used for optimizing smooth, not necessarily convex, 
%   functions, for which the convergence rate is generically superlinear.
%   Surprisingly, it often works very well for functions that are nonsmooth at
%   their minimizers, although the convergence rate is linear at best, and
%   the final inverse Hessian approximation is typically very ill conditioned.
%   Written by Michael Overton, overton@cs.nyu.edu

% parameter defaults
if nargin == 0
   error('bfgs: "pars" is a required input parameter')
end
if nargin == 1
   options = [];
end
%keyboard
options = setdefaults(pars,options);  % fields for all methods
cpufinish = cputime + options.cpumax;
prtlevel = options.prtlevel;
options = setwolfedefaults(options, 0);  % fields for Wolfe line search
if ~isfield(options,'H0')
    options.H0 = []; % not the identity, see bfgs1run.m
end
x0 = options.x0;
nstart = size(x0,2);
for run = 1:nstart
    if prtlevel > 0 & nstart > 1
        fprintf('bfgs: starting point %d\n', run)
    end
    options.cpumax = cpufinish - cputime; % time left
    [x(:,run), f(run), g(:,run), H{run}] = bfgs1run(x0(:,run), pars, options);
    if cputime > cpufinish
        break
    end
end
if nstart == 1
    H = H{1};
end