function ti = cheb_mesh(M,t_lims)
%CHEB_MESH Create Chebyshev maxima grid
% Input:
%   M: number of gridpoints
%   t_lims: lower and upper limit of mesh (default [0 1])
% Output:
%   ti: extrema of T_{M-1} Chebyshev polynomial between t_lims

if nargin < 2
    ti = 1/2*(1-cos(pi/(M-1)*(0:M-1)));
else
    ti = (t_lims(2)-t_lims(1))/2*(1-cos(pi/(M-1)*(0:M-1))) + t_lims(1);
end

end

