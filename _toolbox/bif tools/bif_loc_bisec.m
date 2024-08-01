function [orb1,ds,err] = bif_loc_bisec(orb0,dy0,ds0,funcs,sys,...
    bif_type,opts,bifs)
%BIF_LOC_BISEC Bisection based bifurcation location routine for making
%continuation run endpoints more accurate
% Input:
%   orb0: periodic orbit just befor the bifurcation has been detected
%    -> sig: solution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   dy0: solution tangent vector at orb0
%   ds0: stepsize used to identify a new orbit with a bifurcation
%       condition validated starting from orb0
%   funcs: a function, which returns the governing NAE and its Jacobian wrt
%       state variables and continuation parameters: [F(x), JF(x)]
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
%   bif_type: type of bifurcation expected
%       1) change in stability
%       2) vanishing segment
%       3) grazing
%       4) sliding
%       5) user defined monitor function zero crossing
%   opts: numerical method parameters
%    -> pi: indicies of continuation parameters (length of 1 or 2)
%    -> psa: pseudo-arclength method parameters
%      -> ds_lim: minimum and maximum allowed stepsize [ds_min ds_max]
%           (default [1e-3 10])
%    -> stop: stopping conditions for terminating the continuation run
%       -> pi: indicies of continuation parameters (length of 1 or 2)
%      -> va_tol: tolerance for aborting continuation runs when segment 
%           lengths become too small (default 1e-3)
%      -> gr_tol: tolerance for grazing detection (default 1e-10)
%    -> nr: parameters of the employed Newton Rhapson iteration
%      -> logs: if true print progress of Newton iteration (default false)
%      -> reltol: iteration stopping condition in step norm (default 1e-7)
%      -> abstol: iteration stopping condition in error norm (default 1e-10)
%      -> maxiter: maximum number of iteration steps  (default 10)
%      -> plot: if true plot progress of error function (default false)
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding
%    -> ind: index of bifurcation event in the solution signature
% Output:
%   orb1: periodic orbit at the bifurcation event
%    -> sig: solution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   ds: used pseudo arclength stepsize at orb0
%   err: mp-bvp solution error for orb1

% Initialization
N = length(orb0.sig);   % number of smooth segments
pind = opts.pi;         % indices of continuation parameters
lp = length(pind);      % number of continuation parameters
ds = ds0/2;             % starting stepsize in bisection algorithm
ds_step = abs(ds);      % starting stepsize bisection step
ds_step_min = opts.psa.ds_lim(1)/10; % stopping condition
opts.psa.pi = pind;

% Bisection loop
y0 = [orb0.U; orb0.T; orb0.p(pind).']; % fix starting point
if bif_type == 5
    q0 = feval(sys.q,y0,orb0,sys,pind);
else
    q0 = [];
end
if bif_type == 1 && ~isfield(orb0,'mu_crit')
   [~, orb0.mu_crit, ~] = orb_stab(orb0,sys);
end

while ds_step > ds_step_min

    % Find a new orbit with a psa step
    opts.psa.ds = ds;
    [y1,~,~,err] = psa_step(y0,dy0,funcs,opts);
    orb1 = orb0;
    orb1.U = y1(1:end-N-lp);
    orb1.T = y1(end-N-lp+1:end-lp);
    orb1.p(pind) = y1(end-lp+1:end);

    % Evaluate bifurcation condition for orb1
    switch bif_type
        case 1 % Stability change
            [orb1.mu, orb1.mu_crit, ~] = orb_stab(orb1,sys);
            if (1+opts.nr.abstol-abs(orb0.mu_crit))*...
                (1+opts.nr.abstol-abs(orb1.mu_crit))<0
                bifc = true;
            else
                bifc = false;
            end
        case 2 % Vanishing segment
             if ~isempty(find(orb1.T<opts.stop.va_tol*max(orb1.T),1))
                 bifc = true;
             else
                 bifc = false;
             end
        case 3 % Grazing bifurcation
            if nargin < 8
                [bifc,~,~] = po_graze_det(y0,y1,orb0,sys,opts.stop);
            else
                [bifc,~,~] = po_graze_det(y0,y1,orb0,sys,opts.stop,bifs);
            end
        case 4 % Sliding bifurcation
            if nargin < 8
                [bifc,~] = po_slide_det(y1,orb0,sys,opts.stop);
            else
                [bifc,~] = po_slide_det(y1,orb0,sys,opts.stop,bifs);
            end
        case 5 % User defined monitor function
            orb1.q = feval(sys.q,y1,orb0,sys,pind);
            if any(q0 .* orb1.q < 0)
                bifc = true;
            else
                bifc = false;
            end
    end

    % Update stepsize based on bifc
    ds_step = ds_step/2;
    if bifc
        ds = ds-sign(ds)*ds_step; % reduce if a bifurcation is detected
    else
        ds = ds+sign(ds)*ds_step; % increase if no bifurcation is detected
    end
    
    % disp(ds)
end

end

