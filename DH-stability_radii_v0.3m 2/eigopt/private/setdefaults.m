function options = setdefaults(pars,options)
%  call: options = setdefaults(pars,options)
%  check that fields of pars and options are set correctly and
%  set basic default values for options that are common to various 
%  optimization methods, including bfgs, cgprfr and gradsamp
if ~isfield(pars, 'nvar')
   error('setdefaults: input "pars" must have a field "nvar" (number of variables) is a required parameter')
elseif round(pars.nvar) ~= pars.nvar | pars.nvar <= 0
   error('setdefaults: input "pars.nvar" (number of variables) must be a positive integer')
end
if ~isfield(pars, 'fgname')
   error('setdefaults: input "pars" must have a field "fgname" (name of m-file computing function and gradient)')
end
if ~isfield(options, 'x0')
    options.x0 = [];
end
if isempty(options.x0)
    if isfield(options, 'nrand')
        if round(options.nrand) ~= options.nrand | options.nrand <= 0
            error('setdefaults: input "options.nrand" must be a positive integer when "options.x0" is not provided')
        else
            options.x0 = randn(pars.nvar, options.nrand);
        end
    else
        options.x0 = randn(pars.nvar, 10);
    end
else
    if size(options.x0,1) ~= pars.nvar
        error('setdefaults: input "options.x0" must have "pars.nvar" rows')
    end
    if isfield(options, 'nrand')
        if ~isnonnegint(options.nrand)
            error('setdefaults: input "options.nrand" must be a nonnegative integer')
        else
            options.x0 = [options.x0  randn(pars.nvar, options.nrand)];
        end
    end % no else part, options.x0 is as provided
end
options.nrand = 0; % otherwise MORE random starts may be added by another call
if isfield(options, 'maxit')
    if ~isnonnegint(options.maxit)
        error('setdefaults: input "options.maxit" must be a nonnegative integer')
    end
else
    options.maxit = 100;
end
if isfield(options, 'normtol')
%	keyboard
    if ~isposreal(options.normtol)
        error('setdefaults: input "options.normtol" must be a positive real scalar')
    end
else
    options.normtol = 1.0e-7;
end
if isfield(options, 'fvalquit')
    if ~isreal(options.fvalquit)
        error('setdefaults: input "options.fvalquit" must be a real scalar')
    end

    [nr,nc] = size(options.fvalquit);
    if nr*nc ~= 1
	error('setdefaults: input "options.fvalquit" must be a real scalar')
    end
else
    options.fvalquit = -inf;
end
if ~isfield(options, 'cpumax')
    options.cpumax = inf;
end
if ~isfield(options, 'prtlevel')
    options.prtlevel = 1;
end
