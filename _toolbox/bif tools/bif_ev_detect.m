function [type,orb1] = bif_ev_detect(si,y0,y1,orb0,sys,opts,bifs)
%BIF_EV_DETECT Bifurcation and event detection routines for
% "br12_cont_adapt.m", "br12_cont_fix.m", or "br12_cont_nat.m"
% Input:
%   si: index of the last continuation step
%   y0: previous state vector with its bifurcation parameter [u0; Ti0; pi0]
%   y1: current state vector with its bifurcation parameter [u1; Ti1; pi1]
%   orb0: periodic orbit data structure from the previous orbit
%    -> sig: solution signature (event list)
%    -> U: solution vectors (M*N*n)
%    -> T: solution segment lengths (N)
%    -> p: system parameter vectors (l)
%    -> mu_crit: critical multiplier
%    -> q: auxiliary monitor function evaluations (if q is included in sys)
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> q: monitor function to track during continuation (optional)
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%    -> tau_no: number of distinct time delays
%    -> aut: if true the sytem is truly autonomous, thus the vector field 
%       condition for the final point of the orbit is included 
%       during mondoromy matrix formulation (default true)
%   opts: numerical method parameters
%    -> pi: indicies of continuation parameters (length of 1 or 2)
%    -> psa: pseudo-arclength method parameters
%      -> stab_eval: if true evaluate the stability of all found orbits
%           (default true, can be turned of for reduced calculation times)
%    -> stop: stopping conditions for the continuation run 
%      -> stab_ch: if true stop on changes in stability (default false)
%      -> ext_gr: stop if external grazing is detected (default true)
%      -> int_gr: stop if internal grazing is detected (default true)
%      -> slide: stop if a sliding event is detected (default true)
%      -> user_q: stop on sign changes of the monitor function values 
%           (default false)
%      -> va_tol: tolerance for aborting continuation runs when segment 
%           lengths become too small (default 1e-3)
%      -> gr_tol: tolerance for grazing detection (default 1e-10)
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding
%    -> ind: index of bifurcation event in the solution signature
% Output:
%   type: detected type of bifurcation
%       0) no monitor condition was violated
%       1) change in stability
%       2) vanishing segment
%       3) grazing
%       4) sliding
%       5) user defined monitor function zero crossing
%   orb1: orbit data to be included in the output structure
%    -> U: solution vectors (M*N*n)
%    -> T: solution segment lengths (N)
%    -> p: system parameter vectors (l)
%    -> mu_crit: critical multiplier
%    -> mu: Floquet multipliers
%    -> q: auxiliary monitor function evaluations
%    -> bif_type: bifurcation type if applicable, (possible types: 
%       stability change, grazing, sliding, vanishing segment, user 
%       defined zero crossing)

% Initialization
N = length(orb0.sig);    % number of segments
lp = length(opts.pi);    % number of continuation parameters
type = 0;                % solution bifurcation flag (empty by default)

% Fill up a temporary orbit data structure
orb1 = orb0;
orb1.U = y1(1:end-N-lp);
orb1.T = y1(end-N-lp+1:end-lp); 
orb1.p(opts.pi) = y1(end-lp+1:end);
orb1.bif_type = [];

% Evaluate orbit stability
if opts.psa.stab_eval
    [orb1.mu, orb1.mu_crit, ~] = orb_stab(orb1,sys);
end

% Auxiliary monitor functions
if isfield(sys,'q')
    orb1.q = feval(sys.q,y1,orb0,sys,opts.pi);
    if any(orb1.q .* orb0.q<0)
        text = 'Sign change in user defined monitor function';
        if ~isempty(orb1.bif_type)
            orb1.bif_type = [orb1.bif_type '; ' text];
        else
            orb1.bif_type = text;
        end
        fprintf('   -> User defined zero crossing detected at step %i\n',si);
        if opts.stop.user_q
            type = [type 5];
        end
    end
else
    orb1.q = [];
end

% Check for sliding events
if nargin < 7
    [slide,s_ind] = po_slide_det(y1,orb0,sys,opts.stop);
else
    [slide,s_ind] = po_slide_det(y1,orb0,sys,opts.stop,bifs);
end
if slide
    text = sprintf('Boundary sliding (at event %i)',s_ind);
    if ~isempty(orb1.bif_type)
        orb1.bif_type = [orb1.bif_type '; ' text];
    else
        orb1.bif_type = text;
    end
    fprintf('   -> Sliding bifurcation detected at step %i\n',si);
    if opts.stop.slide
        type = [type 4];
    end
end

% Check for grazing events
if nargin < 7
    [graze,g_int,g_ind] = po_graze_det(y0,y1,orb0,sys,opts.stop);
else
    [graze,g_int,g_ind] = po_graze_det(y0,y1,orb0,sys,opts.stop,bifs);
end
if graze
    if g_int
        text = sprintf('Interior grazing (in event condition %i)',g_ind);
    else
        text = sprintf('Exterior grazing (at event %i)',g_ind);
    end
    if ~isempty(orb1.bif_type)
        orb1.bif_type = [orb1.bif_type '; ' text];
    else
        orb1.bif_type = text;
    end
    fprintf('   -> Grazing bifurcation detected at step %i\n',si);
    if (opts.stop.int_gr && g_int) || (opts.stop.ext_gr && ~g_int)
        type = [type 3];
    end
end

% Check orbits for vanishing segments
vi = find(orb1.T<opts.stop.va_tol*max(orb1.T)); % index of vanishing segments
if ~isempty(vi)
    text = sprintf('Vanishing segment at index %i',vi(1));
    if ~isempty(orb1.bif_type)
        orb1.bif_type = [orb1.bif_type '; ' text];
    else
        orb1.bif_type = text;
    end
   fprintf('   -> Vanishing segment detected at step %i\n',si);
   type = [type 2];
end

% Mark changes in stability
if opts.psa.stab_eval
    if (1+opts.nr.abstol-abs(orb0.mu_crit))*...
        (1+opts.nr.abstol-abs(orb1.mu_crit))<0
        if abs(imag(orb1.mu_crit))>opts.nr.abstol
            text = 'Stability change (Hopf)';
        elseif real(orb1.mu_crit)>0
            text = 'Stability change (Saddle node)';
        else
            text = 'Stability change (Period doubling)';
        end
        if ~isempty(orb1.bif_type)
            orb1.bif_type = [orb1.bif_type '; ' text];
        else
            orb1.bif_type = text;
        end
        fprintf('   -> Change in stability detected at step %i\n',si);
        % Stop continuation if necessary
        if opts.stop.stab_ch
            type = [type 1];
        end
    end
end

% Preoritize between events
% user defined monitor function > sliding > grazing > vanishing > stability
type = max(type);

end

