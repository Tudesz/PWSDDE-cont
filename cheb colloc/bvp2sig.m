function [t,u] = bvp2sig(U,T,M,res)
%BVP2SIG Convert solution data from BVP vector to form signal form
% Input:
%   U: solution in column vector form
%   T: lengths of smooth segments (N)
%   M: interpolation resolution
%   res: output resolution (default M -> no interpolation)
% Output:
%   t: time mesh of solution (M*N)
%   u: signal on t (n,M*N)

% Initialization
N = length(T);
n = length(U)/(N*M);

if nargin<4
% Generate output
    t = zeros(1,M*N);
    u = zeros(n,M*N);
    for i = 1:N
        ii = (i-1)*M+1:i*M; % indicies of u and t
        ii0 = (i-1)*M*n+1:i*M*n; % indicies of u0
        u(:,ii) = reshape(U(ii0),M,n).';
        t(ii) = cheb_mesh(M,[sum(T(1:i-1)) sum(T(1:i))]);
    end
else
    t = zeros(1,res*N);
    u = zeros(n,res*N);
    for i = 1:N
        ii = (i-1)*res+1:i*res;
        t(ii) = linspace(sum(T(1:i-1)),sum(T(1:i)),res);
        for j = 1:res
            [u(:,ii(j)),~,~] = po_interp(t(ii(j)),U,T,M);
        end
    end


end

