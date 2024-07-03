function [uq,iq,fiq] = po_interp(tq,u0,Ti,M)
%PO_INTERP Interpolation on the multi segment Chebyshev mesh at tq
% Input:
%   tq: querry point (0<tq<T)
%   u0: solution in column vector form
%   Ti: lengths of smooth segments (N)
%   M: interpolation resolution
% Output:
%   uq: state at tq (n)
%   iq: segment index
%   fiq: lagrange interpolation coefficients for tq

% Initialization
N = length(Ti);
n = length(u0)/(N*M);

% Find correspoding solution segment
% if tq < 0 || tq > sum(Ti)
%     warning('Request outside orbit boundaries at tq = %0.3f',tq)
% end
iq = find(cumsum(Ti)>=tq,1);

% Interpolate state on segment iq
tq_iq = (tq-sum(Ti(1:iq-1)))/Ti(iq); % rescaled querry point (in [0 1])
fiq = lag_coeff(M,tq_iq); % lagrance coefficients on [0 1] Chebishev mesh
u0_iq = u0((iq-1)*M*n+1:iq*M*n); % state vectors on segment iq
uq = (fiq*reshape(u0_iq,M,n)).';

end

