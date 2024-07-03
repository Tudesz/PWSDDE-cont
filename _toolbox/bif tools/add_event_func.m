function X = add_event_func(U,T,p,orb,sys,bif_ind,lp,type,ev_id)
%ADD_EVENT_FUNC Additional event condition for the continuation of double
%events in PWS-DDE periodic orbits
% Input:
%   U: free state variables (N*M*n)
%   T: free segment lengths (N)
%   p: parameter vector
%   orb: periodic orbit data structure (only the metadata part is used)
%    -> sig: sloution signature (event list)
%    -> n: number of degrees of freedom
%    -> M: Chebyshev mesh resolution
%   sys: names of the functions that define the system
%    -> f: vector field and its Jacobians
%    -> e: event function, map and corresponding Jacobians
%    -> tau: time delay and its parameter Jacobian
%    -> tau_no: number of distinct time delays
%   bif_ind: index of the bifurcation in the solution signature
%   lp: bifurcation parameter index (default all)
%   type: flag for which property is requested:
%       1, The function evaluation itself
%       2, Jacobian wrt U and T
%       3, Parameter Jacobian
%   ev_id: index of the extra event condition in sys.e
% Output:
%   X: appropriate evaluation of the event functions of Jacobians

% Initialization
M = orb.M;                  % mesh resolution
n = orb.n;                  % system dimension
N = length(orb.sig);        % number of segments
D0 = cheb_diff(M);          % base differentiation matrix
[~,x0] = bvp2sig(U,T,M);    % unpack solution vector
x = x0(:,bif_ind*M);        % state at sliding event

% Interpolation of delayed terms
[x0_tau,i_tau,fi_tau,~,dT] = po_delay_interp(U,T,p,M,sys);

% Event condition and delay function initialization
xd = squeeze(x0_tau(:,:,bif_ind*M));
ev_cond = @(t,l) feval(sys.e,x,xd,p,ev_id,t,l);
tau_f = @(k,t) feval(sys.tau,p,k,t);

% Return requested output
switch type
    
    case 1
        % Event condition
        X = ev_cond(1,0);

    case 2
        % Event condition Jacobain
        X = zeros(1,n*N*M+N);
        X(1,(bif_ind-1)*M*n+M:M:bif_ind*M*n) = ev_cond(2,0);
        % Account for delayed terms
        for k = 1:sys.tau_no
            % Interpolation
            fi_k = squeeze(fi_tau(bif_ind*M,k,:)).';
            dfi_k = fi_k*D0; 
            i_k = i_tau(bif_ind*M,k);
            ik = (i_k-1)*M*n+1:i_k*M*n;
            dT_k = squeeze(sum(dT(bif_ind*M,k,:,:),4)).';
            % Place elements in Jacobian matrix
            Jev_k = ev_cond(2,k);
            X(1,ik) = X(1,ik) + kron(Jev_k,fi_k);
            X(1,end-N+1:end) = X(1,end-N+1:end) ...
                + kron(Jev_k,dfi_k)*U(ik)*dT_k;
            if tau_f(k,1) == 2 % extra term for neutral delays
                X(1,end-N+i_k) = X(1,end-N+i_k) - ...
                    kron(Jev_k,fi_k)*U(ik)/T(i_k);
            end
        end

    case 3
        % Event condition parameter Jacobian
        X = ev_cond(3,0);
        % Account for delayed terms
        for k = 1:sys.tau_no
            fi_k = squeeze(fi_tau(sl_ind*M,k,:)).';
            dfi_k = (D0.'*fi_k.').';
            i_k = (i_tau(sl_ind*M,k)-1)*M*n+1:i_tau(sl_ind*M,k)*M*n;
            dtau_k = -1/T(i_tau(sl_ind*M,k))*tau_f(k,3);
            X = X + kron(ev_cond(2,k),dfi_k)*U(i_k)*dtau_k;
        end
        % Truncate if necessary
        if ~isempty(lp)
            X = X(lp);
        end
end

end

