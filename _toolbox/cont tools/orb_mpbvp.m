function [F,JF] = orb_mpbvp(x,orb,sys,pind,bifs)
%ORB_MPBVP Evaluate the NAE and its Jacobian corresponding to the MP-BVP
%of a periodic orbit
%   orb: starting periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%   pind: indices of continuation parameters (can be empty)
%   bifs: extra data for finding bifurcation points (optional input)
%    -> type: 1) grazing 2) sliding 3) user defined bifurcation
%    -> ind: index of the bifurcation event in the solution signature
%    -> pi: index of a free system parameter to be corrected
%    -> f: function describing the user defined zero condition and 
%       its Jacobains (only required if bifs.type==3)
% Output:
%   func: function returning the governing NAE and its appropriate Jacobian
%       f(x,p) = [F(x,p), JF(x,p)]

if nargin < 4
    pind = []; % no continuation variable by default
end

% Extract necessary orbit data
N = length(orb.sig); % number of smooth segments
U = x(1:end-N-length(pind)); % state vector
T = x(end-N+1-length(pind):end-length(pind)); % segment lengts
par = orb.p; % parameter vector (default)
par(pind) = x(end-length(pind)+1:end); % insert the varied parameters in their place

% Interpolate all delayed terms
del = po_delay_interp(U,T,par,orb.M,sys); 

if isempty(pind)
    % Correction of a single regular point
    F = mpbvp(U,T,par,orb,sys,del);
    JF = mpbvp_Ju(U,T,par,orb,sys,del);

elseif length(pind) == 1 && nargin < 5
    % Regular point on a one parameter branch
    F = mpbvp(U,T,par,orb,sys,del);
    JU = mpbvp_Ju(U,T,par,orb,sys,del);
    JP = mpbvp_Jp(U,T,par,orb,sys,del);
    JF = [JU JP(:,pind)];

elseif length(pind) <= 2
    % Bifurcation points correction/continuation
        switch bifs.type
            case 1 % Grazing bifurcation point
                F = [mpbvp(U,T,par,orb,sys,del);...
                    mpbvp_gr(U,T,par,orb,sys,del,bifs.ind)];
                JU = [mpbvp_Ju(U,T,par,orb,sys,del); ...
                    mpbvp_gr_Ju(U,T,par,orb,sys,del,bifs.ind)];
                JP = [mpbvp_Jp(U,T,par,orb,sys,del); ...
                    mpbvp_gr_Jp(U,T,par,orb,sys,del,bifs.ind)];
            case 2 % Sliding bifurcation point
                 F = [mpbvp(U,T,par,orb,sys,del);...
                    mpbvp_sl(U,T,par,orb,sys,del,bifs.ind)];
                JU = [mpbvp_Ju(U,T,par,orb,sys,del); ...
                    mpbvp_sl_Ju(U,T,par,orb,sys,del,bifs.ind)];
                JP = [mpbvp_Jp(U,T,par,orb,sys,del); ...
                    mpbvp_sl_Jp(U,T,par,orb,sys,del,bifs.ind)];
            case 3 % User defined bifurcation condition
                F = [mpbvp(U,T,par,orb,sys,del);...
                    bifs.f(U,T,par,orb,sys,del,bifs.ind,1)];
                JU = [mpbvp_Ju(U,T,par,orb,sys,del); ...
                    bifs.f(U,T,par,orb,sys,del,bifs.ind,2)];
                JP = [mpbvp_Jp(U,T,par,orb,sys,del); ...
                    bifs.f(U,T,par,orb,sys,del,bifs.ind,3)];
        end
    JF = [JU JP(:,pind)];

else
   error('Only one and two parameter continuation is supported');
end


end

