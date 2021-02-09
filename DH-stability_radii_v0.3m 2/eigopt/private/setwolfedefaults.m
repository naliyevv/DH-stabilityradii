function options = setwolfedefaults(options, CG)
%  check Wolfe line search fields for options and set defaults
%  CG is 1 for CG methods, 0 for BFGS methods
if isfield(options, 'strongwolfe')
    if options.strongwolfe ~= 0 & options.strongwolfe ~= 1
        error('setwolfedefaults: input "options.strongwolfe" must be 0 or 1')
    end
else
    if CG
        options.strongwolfe = 1;  % needed for convergence analysis
    else
        options.strongwolfe = 0;  % not needed for convergence analysis
        % strong Wolfe is very complicated and is bad for nonsmooth functions
    end
end
if isfield(options, 'wolfe1') 
    if ~isposreal(options.wolfe1)
        error('setwolfedefaults: input "options.wolfe1" must be a positive real scalar')
    end
else
    options.wolfe1 = 1e-4;
end
if isfield(options, 'wolfe2')
    if ~isposreal(options.wolfe2)
        error('setwolfedefaults: input "options.wolfe2" must be a positive real scalar')
    end
elseif CG == 1
    options.wolfe2 = 0.49;  % must be < .5 for CG convergence theory
else
    options.wolfe2 = 0.9;   % must be < 1 for BFGS update to be pos def
end
if options.wolfe1 <= 0 | options.wolfe1 >= options.wolfe2 | options.wolfe2 >= 1
    fprintf('setwolfedefaults: Wolfe line search parameters violate usual requirements')
end
if options.prtlevel > 0
    if options.strongwolfe & ~CG
        fprintf('Strong Wolfe line search selected, but for BFGS or LMBFGS\n')
        fprintf('weak Wolfe may be preferable, especially if f is nonsmooth\n')
    elseif ~options.strongwolfe & CG 
        fprintf('Weak Wolfe line search selected, but for CG\n')
        fprintf('this often fails: use strong Wolfe instead')
    end
end
            