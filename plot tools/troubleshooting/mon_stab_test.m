function [orb_p0, orb_pT] = mon_stab_test(orb,sys,prt)
%MON_STAB_TEST Test how the approximate monodromy matrix maps perturbations
% forward in time
% Input:
%   orb: periodic orbit data structure
%    -> U: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p: parameter vector (l)
%    -> sig: solution signature (n)
%   sys: names of the functions that define the 
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> mode_no: number of distinct vector field modes
%    -> event_no: number of disctinc events
%   prt: norm of applied perturbation or a perturbation vector of size M*n
%       or M*n+N (default 0)
% Output:
%   orb_p0: perturbed orbit data structure
%   orb_pT: new orbit with perturbations mapped forward in time using the 
%       approximate monodromy matrix

% Initialize orbit parameters
U = orb.U;
T = orb.T;

% apply some perturbation
if length(prt) == 1
    yp = [prt*norm(U)*randn(size(U)); zeros(size(T))];
elseif length(prt) == length(U)
    yp = [prt; zeros(size(T))];
elseif length(prt) == length(U) + length(T)
    yp = prt;
end
orb_p0 = orb;
orb_p0.U = U + yp(1:length(U));
orb_p0.T = T + yp(length(U)+1:end);

% forward propagation via the monodromy matrix
[mu,~,Th,~] = orb_stab(orb,sys);
K = length(yp); 
ntau = length(mu)/K;
ypT = Th(end-K+1:end,:)*repmat(yp,ntau,1);
orb_pT = orb;
orb_pT.U = U + ypT(1:length(U));
orb_pT.T = T + ypT(length(U)+1:end);


end

