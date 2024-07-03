function F = cheb_int(ti,fi)
%CHEB_INT Calculate Clenshaw-Curtis quadrature for integrating a function 
% defined on the Chebyshev mesh
% Input:
%   ti: Chebyshev mesh where the function values are known
%   fi: function values at ti
% Output:
%   F: definite intergal of f on ti

T = ti(end)-ti(1);
wi = cheb_quad(length(ti));
F = (T*wi*fi).';

end

