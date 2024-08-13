function [type,conv,bif_type] = br12_term_err(si,y1,err1,orb0,opts)
%BR12_TERM_ERR Terminal error/event handling for continuation runs done via
% "br12_cont_adapt.m", "br12_cont_fix.m", or "br12_cont_nat.m"
% Input:
%   si: index of the last continuation step
%   y1: current state vector with its bifurcation parameter [u1; Ti1; pi1]
%   err1: error of the last pseudo arclength step
%   orb0: periodic orbit data structure from the previous orbit
%    -> sig: solution signature (event list)
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> q: monitor function to track during continuation (optional)
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%   opts: numerical method parameters
%    -> pi: indicies of continuation parameters (length of 1 or 2)
%    -> stop: stopping conditions for the continuation run 
%      -> conv: if true stop when Newton iteration fails to converge
%           (default true)
%      -> Tneg: stop when a negative segment lenght is encountered 
%           (default true)
% Output:
%   type: detected type of bifurcation
%      -1) terminal events with the last point saved
%      -2) terminal events with the last point omitted
%   conv: flag for marking solver convergence
%   bif_type: terminal event type, if applicable (possible types: 
%       negative segment, parameter boundary)

% Initialization
N = length(orb0.sig);    % number of segments
lp = length(opts.pi);    % number of continuation parameters
type = 0;                % solution bifurcation flag (empty by default)
bif_type = [];           % no output message by default
  
% Warn/stop in case of non-convergent solutions
conv = norm(err1) < opts.nr.abstol*1e3;
if ~conv
    warning('Solution in br12_cont_adapt did not converge at step %i',si);
    bif_type = sprintf('Non-convergent solution!');
    if opts.stop.conv
        type = -2;
    end
end

% Warn/stop in case of negative segment lengths
Ti = y1(end-N-lp+1:end-lp);
if type > -2 && (any(Ti<0,'all') || any(isnan(Ti),'all'))
    warning('Negative segment length encountered at step %i',si);
    [~,vi] = min(Ti);
    bif_type = sprintf('Negative segment at index %i',vi);
    if opts.stop.Tneg
        type = -1;
    end
end

end

