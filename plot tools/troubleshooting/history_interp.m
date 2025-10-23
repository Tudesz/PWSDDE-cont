function uq = history_interp(orb,tq)
%HISTORY_INTERP  Interpolation on the multi segment Chebyshev mesh at tq
%for the evaluation of initial history data
% Input:
%   orb: periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   tq: querry points on [-tau 0]
% Output:
%   uq: state at tq (n)

nt = length(tq);
uq = zeros(orb.n,nt);
T = sum(orb.T);
for i = 1:nt
    uq(:,i) = po_interp(ceil(-tq(i)/T)*T+tq(i),orb.U,orb.T,orb.M);
end

end

