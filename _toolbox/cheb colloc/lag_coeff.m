function li = lag_coeff(M,tq)
%LAG_COEFF2 Coefficients of Lagrange interpolation on the [0 1]
%Chebyshev type 2 grid using the barycentric form
% Input:
%   M: number of mesh points
%   tq: query point
% Output:
%   li: Lagrange interpolation coefficient li(tq)

% Barycentric weights
wk = (-1).^(1:M);
wk(1) = wk(1)/2; 
wk(end) = wk(end)/2;

% Chebyshev mesh
tk = cheb_mesh(M);

% Interpolation coefficients
if any(tk==tq)
    li = zeros(1,M);
    li(tk==tq) = 1;
else
    li = wk./(tq-tk);
    li = li./sum(li);
end

end

