function [alpha, xalpha, falpha, gradbeta, fail, beta, nfeval] = ...
          linesch_ww(x0, f0, g0, d, pars, c1, c2, prtlevel);
% Copyright Michael Overton
% LINESCH_WW Line search enforcing weak Wolfe conditions, suitable
%            for minimizing both smooth and nonsmooth functions
% call:  [alpha, xalpha, falpha, gradbeta, fail, beta, nfeval] = ...
%         linesch_ww(x0, f0, g0, d, pars, c1, c2, prtlevel);
%  Input
%   x0:      intial point
%   f0:      function value at x0
%   g0:      directional derivative of function at x0 in direction d
%   d:       search direction  
%   pars:    a structure that specifies the function name as well
%            anything else that the user needs to access in programming the
%            function and gradient values
%        pars.fgname:  name of function that returns function and gradient
%            it expects as input only x and pars, a parameter structure 
%            it is invoked by: [f,g] = feval(fgname, x, pars)
%   c1: Wolfe parameter for the sufficient decrease condition 
%          f(x0 + t d) <= f0 + c1*t*g0   (default 1e-4)
%   c2: Wolfe parameter for the WEAK condition on directional derivative
%          (grad f)(x0 + t d)'*d >= -c2*g0   (default 0.9)
%        these should normally satisfy 0 < c1 < c2 < 1, but may
%        want c1 = c2 = 0 when function known to be nonsmooth
%        note: setting c2 very small may interfere with superlinear
%           convergence when function is smooth
%   prtlevel: 0 for no printing, 1 minimal (default), 2 verbose 
%
%  Output:
%   alpha:   steplength satisfying weak Wolfe conditions if one was found,
%             otherwise left end point of interval bracketing such a point
%             (possibly 0)
%   xalpha:  x0 + alpha*d
%   falpha:  f(x0 + alpha d)
%   gradbeta:(grad f)(x0 + beta d)  (see below for definition of beta)
%             (nan if beta = inf)
%   fail:    0 if both Wolfe conditions satisfied
%            1 if one or both Wolfe conditions not satisfied but an
%               interval was found bracketing a point where both satisfied
%           -1 if no such interval was found, f may be unbounded below
%   beta:    steplength satisfying weak Wolfe conditions if one was found,
%             otherwise right end point of interval bracketing such a point
%             (inf if no such finite interval found)
%   nfeval:  number of function evaluations
%  Note: the bracketing property assumes f is continuous, but not that
%  it is smooth (though we probably need to assume smooth a.e.).
%  The gradient at alpha, not beta, is returned, because in the case that
%  alpha < beta, it is normally the gradient at beta that is needed. 
%  The strange outpout arg order is to be consistent with linesch_sw. 

% The weak Wolfe line search is far less complicated that the standard 
% strong Wolfe line search that is discussed in many texts. It appears
% to have no disadvantages compared to strong Wolfe when used with
% Newton or BFGS methods on smooth functions, and it is essential for 
% the application of BFGS or bundle to nonsmooth functions as done in HANSO.
% However, it is NOT recommended for use with conjugate gradient methods,
% which require a strong Wolfe line search for convergence guarantees.
% Weak Wolfe requires two conditions to be satisfied: sufficient decrease
% in the objective, and sufficient increase in the directional derivative
% (not reduction in its absolute value, as required by strong Wolfe).
%
% There are some subtleties for nonsmooth functions.  In the typical case
% that the directional derivative changes sign somewhere along d, it is
% no problem to satisfy the 2nd condition, but descent may not be possible
% if the change of sign takes place even when the step is tiny. In this
% case it is important to return the gradient corresponding to the positive 
% directional derivative even though descent was not obtained. On the other 
% hand, for some nonsmooth functions the function decrease is steady
% along the line until at some point it jumps to infinity, because an
% implicit constraint is violated.  In this case, the first condition is
% satisfied but the second is not. All cases are covered by returning
% the end points of an interval [alpha, beta] which is constructed to
% contain a point satisfying the weak Wolfe conditions assuming f is 
% continuous (and smooth almost a.e.?), and returning the function value
% at alpha, but the gradient at beta.
%
% For functions that are known to be nonsmooth, setting the second Wolfe
% parameter to zero makes sense, especially for a bundle method.  However,
% for smooth functions this may prevent superlinear convergence.
%
%  Written by Michael Overton (overton@cs.nyu.edu)
if nargin < 6  % check if the optional Wolfe parameters were passed
    c1 = 1e-4;
end
if nargin < 7
    c2 = 0.9;
end
if nargin < 8
    prtlevel = 1;
end
if c1 < 0 | c1 > c2 | c2 >= 1 % allows c1 = c2 = 0
   if prtlevel > 0
       fprintf('linesch_ww: Wolfe parameters do not satisfy 0 <= c1 <= c2 <= 1')
   end
end
fgname = pars.fgname;
alpha = 0;  % lower bound on steplength conditions
xalpha = x0;
falpha = f0;
beta = inf;  % upper bound on steplength satisfying weak Wolfe conditions
gradbeta = nan*ones(size(x0));
if g0 >= 0
    error('linesch_ww: g0 is negative, indicating d not a descent direction')
end
dnorm = norm(d);
if dnorm == 0
    error('linesch_ww: d is zero')
end
t = 1;  % important to try steplength one first
nfeval = 0;
nbisect = 0;
nexpand = 0;
% the following limits are rather arbitrary
% nbisectmax = 30; % 50 is TOO BIG, because of rounding errors
nbisectmax = max(30, round(log2(1e5*dnorm))); % allows more if ||d|| big
nexpandmax = max(10, round(log2(1e5/dnorm))); % allows more if ||d|| small
done = 0;
while ~done
    x = x0 + t*d;
    nfeval = nfeval + 1;
%	keyboard
    [f,grad] = feval(fgname, x, pars);
    gtd = grad'*d;
    % the first condition must be checked first
    if f > f0 + c1*t*g0 % first condition violated, gone too far
        beta = t;
        gradbeta = grad;
    elseif gtd < c2*g0 % second condition violated, not gone far enough
        alpha = t;
        xalpha = x;
        falpha = f;
    else   % quit, both conditions are satisfied
        fail = 0;
        alpha = t;
        xalpha = x;
        falpha = f;
        beta = t;
        gradbeta = grad;
        return
    end
    % setup next function evaluation
    if beta < inf
        if nbisect < nbisectmax
            nbisect = nbisect + 1;
            t = (alpha + beta)/2; % bisection
        else
            done = 1;
        end
    else
        if nexpand < nexpandmax
            nexpand = nexpand + 1;
            t = 2*alpha;  % still in expansion mode
        else
            done = 1;
        end
    end 
end % loop
% Wolfe conditions not satisfied: there are two cases
if beta == inf % minimizer never bracketed
    fail = -1;
    if prtlevel > 1
        fprintf('Line search failed to bracket point satisfying weak ');
        fprintf('Wolfe conditions, function may be unbounded below\n')
    end
else % point satisfying Wolfe conditions was bracketed
    fail = 1;
    if prtlevel > 1
        fprintf('Line search failed to satisfy weak Wolfe conditions')
        fprintf(' although point satisfying conditions was bracketed\n')
    end
end