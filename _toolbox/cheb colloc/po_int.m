function I = po_int(f,u0,Ti,M)
%PO_INT Integrate a function f over a periodic orbit defined on a 
%piecewise-Chebyshev mesh 
% Input:
%   f: function to integrat f(x) R^n->R^m
%   u0: solution in column vector form
%   Ti: lengths of smooth segments (N)
%   M: interpolation resolution

% Evaluate function values
[ts,us] = bvp2sig(u0,Ti,M);
fs = f(us);

% Integrate segment by segment
I = zeros(size(fs,1),1);
for i = 1:length(Ti)
    ii = (i-1)*M+1:i*M; % indicies of u and t
    I = I + cheb_int(ts(ii),fs(:,ii).');
end

