function [sl,ind] = po_slide_det(y0,orb,sys,opts,bifs)
%PO_SLIDE_DET Sliding event detection in continuation of periodic orbits based
%on reaching +-1 in the sliding surface condition
% Input:
%   y0: previous state vector with its bifurcation parameter [u0; Ti0; pi0]
%   orb: periodic orbit data structure (only the metadata part is used)
%    -> sig: sloution signature (event list)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%   opts: solver options
%    -> pi: index of continuation parameter
%   bifs: extra data for two parameter continuation (optional input)
%    -> type: 1) grazing 2) sliding
%    -> ind: index of bifurcation event in the solution signature
% Output:
%   sl: bolean flag for sliding bifurcations
%   ind: index of event which underwent a sliding event

% Initialization
sl = false;             % no bifurcation by default
ind = [];               % no index by default
M = orb.M;              % mesh resolution
N = length(orb.sig);    % number of segments

% Unpack solution vector
lp = length(opts.pi);
p0 = orb.p; p0(opts.pi) = y0(end-lp+1:end);
Ti0 = y0(end-N-lp+1:end-lp);

% Function definitions
pi_ej = @(j,i) feval(sys.e,[],[],[],orb.sig(j),7,i);    % incoming and outgoing modes at events
f_mj = @(mj,x,xd) feval(sys.f,x,xd,p0,mj,1,0);    % vector field in mode mj
Jh_ej = @(j,x,xd) feval(sys.e,x,xd,p0,orb.sig(j),2,0);  % event condition jacobian at ej

% Evaluate current and delayed terms
[~,x0] = bvp2sig(y0(1:end-N-lp),Ti0,M);
del0 = po_delay_interp(y0(1:end-N-lp),Ti0,p0,M,sys);
x0_tau = del0.ud;

% Span all segment boundaries
for j = 1:N
    % Extract event data
    m_in = pi_ej(j,1); % mode before event
    m_out = pi_ej(j,2); % mode after event
    f_in = f_mj(m_in,x0(:,j*M),squeeze(x0_tau(:,:,j*M))); % vector field before event
    f_out = f_mj(m_out,x0(:,j*M),squeeze(x0_tau(:,:,j*M))); % vector field after event
    dh = Jh_ej(j,x0(:,j*M),squeeze(x0_tau(:,:,j*M))); % event surface jacobian
    
    % Check sliding condition
    if m_in ~= m_out
        h_sl = (dh*f_in + dh*f_out)/(dh*f_in - dh*f_out);
        if abs(h_sl) < 1 && (lp == 1 || bifs.ind ~= j)
            % Omit false sliding events when following bifurcations
            sl = true;
            ind = j;
            break
        end
    end
    
end


end

