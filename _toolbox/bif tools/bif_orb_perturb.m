function orb_p = bif_orb_perturb(orb,sys,per_opts)
%BIF_ORB_PERTURB Perturb an orbit in the appropriate direction for branch
%swithing
% Input:
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
%   per_opts: data on the previously detected bifurcation
%    -> type: 1) period doubling 2) saddle node
%    -> dU: perturbation level in U (default 1)
%    -> dT: perturbation level in T (default 1)
% Output:
%   orb_p: new, perturbed periodic orbit data structure
%    -> sig: sloution signature (event list)
%    -> U: state variable vector (M*N*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution

% Initialization
N = length(orb.sig);
if ~isfield(per_opts,'dU')
    per_opts.dU = 1;
end
if ~isfield(per_opts,'dT')
    per_opts.dT = 1;
end

% Select perturbation type
switch per_opts.type

    case 1 
        % Period doubling bifurcation
        [~,~,~,Vc] = orb_stab(orb,sys); % find the critical eigenvector of the monodromy matrix
        vU = per_opts.dU*Vc(1:end-N); % perturvation in U
        vT = per_opts.dT*Vc(end-N+1:end); % perturbation in T
        orb_p = orb_convert(orb,sys,-1); % double the orbit
        % perturb the first part
        orb_p.T(1:N) = orb.T + vT; 
        orb_p.U(1:length(vU)) = orb.U + vU;
        % perturb the second part
        orb_p.T(N+1:end) = orb.T - vT; 
        orb_p.U(length(vU)+1:end) = orb.U - vU;

    case 2
        % Saddle node bifurcation (transcritical or pitchfork scenario)
        [~,~,~,Vc] = orb_stab(orb,sys); % find the critical eigenvector of the monodromy matrix
        vU = per_opts.dU*Vc(1:end-N); % perturvation in U
        vT = per_opts.dT*Vc(end-N+1:end); % perturbation in T
        orb_p = orb;
        orb_p.T = orb.T + vT; 
        orb_p.U = orb.U + vU;

end

end

