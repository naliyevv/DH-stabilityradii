function inni = isnonnegint(x)
% return true if x is a nonnegative integer, false otherwise 
[nr,nc] = size(x);
if nr*nc ~= 1
    inni = 0;
else % following is OK since x is scalar
    inni = (isreal(x) & round(x) == x & x >= 0);
end
