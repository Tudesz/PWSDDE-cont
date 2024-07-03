function par_1 = pl_insert(par_0,p,ind)
%PL_INSERT Parameters p into indicies ind in par_0
% Input:
%   par_0: default parameter vector
%   p: parameters to insert
%   ind: indicies of parameters (same length as p)
% Output:
%   par_1: updated parameter vector

par_1 = par_0;
par_1(ind) = p;
end

