function wi = cheb_quad(M)
%CHEB_QUAD Clenshaw-Curtis quadrature weights for a Chebyshev grid 
% F(t_M) = sum(wi*f(ti))
% Input:
%   M: number of gridpoints
% Output:
%   wi: quadrature weights to a Chebyshev grid defined on [0 1]

% Clenshaw-Curtis quadrature on [1 -1] 
theta = pi*(1:M-2)/(M-1);
v = ones(1,M);
if mod(M-1,2)==0
    for k=1:(M-1)/2-1
        v(2:M-1) = v(2:M-1) - 2*cos(2*k*theta)/(4*k^2-1); 
    end
    v(2:M-1) = v(2:M-1) - cos((M-1)*theta)/((M-1)^2-1);
    v(1) = 1/((M-1)^2-1);
    v(end) = v(1);
else
    for k=1:(M-2)/2
        v(2:M-1) = v(2:M-1) - 2*cos(2*k*theta)/(4*k^2-1);
    end
    v(1) = 1/((M-1)^2);
    v(end) = v(1);
end
v(2:end-1) = 2*v(2:end-1)/(M-1);

% rescale quadrature weights to the interval [0 1]
wi = v/2;

end
