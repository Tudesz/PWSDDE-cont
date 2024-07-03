function Dm = cheb_diff(M)
%CHEB_DIFF Create Chebyshev spectral differentiation matrix on [0 1]
% Input:
%   M: mesh resolution (number of grid points)
% Output:
%   Dm: spectral differntiation matrix

% Initialization
D = zeros(M);
m = M-1; %indexing goes from 0 to M-1
x = cos(pi/m*(0:m)); % chebyshev mesh on [-1,1]

% Fill up differentiation matrix
D(1,1) = (2*m^2+1)/6;
D(m+1,1) = -1/2*(-1)^m;
D(1,m+1) = -D(m+1,1);
D(m+1,m+1) = -D(1,1);
for i = 2:m
    D(i,1) = -1/2*(-1)^(i-1)/(1-x(i));
    D(1,i) = 2*(-1)^(i-1)/(1-x(i));
    D(i,m+1) = 1/2 *(-1)^(m+i-1)/(1+x(i));
    D(m+1,i) = -2*(-1)^(m+i-1)/(1+x(i));
end
for i = 2:m
    for j = 2:m
        if i == j
            D(i, j) = -x(j)/(2*(1-x(j)^2));
        else
            D(i, j) = (-1)^(i+j-2)/(x(i)-x(j));
        end
    end
end

% Rescale differntiation matrix to [0 1]
Dm = -2*D;

end

