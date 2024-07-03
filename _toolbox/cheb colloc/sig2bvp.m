function [U,T] = sig2bvp(t,u,M)
%SIG2BVP Convert solution data from signal form to BVP vector form
% Input:
%   t: time mesh of solution (M*N)
%   u: signal on t (n,M*N)
%   M: interpolation resolution
% Output:
%   U: solution in column vector form
%   T: lengths of smooth segments (N)

% Initialization
N = length(t)/M;
n = size(u,1);
T = zeros(N,1);
U = zeros(N*M*n,1);

% Generate output
for i = 1:N
    ii = (i-1)*M+1:i*M; % indicies of u and t
    ii0 = (i-1)*M*n+1:i*M*n; % indicies of u0
    U(ii0) = reshape(u(:,ii).',[],1);
    T(i) = t(ii(end))-t(ii(1));
end

end

