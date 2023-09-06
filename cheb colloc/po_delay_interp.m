function [ud,id,fi,nd,dT] = po_delay_interp(u0,Ti,M,tau,res)
%PO_DELAY_INTERP Evaluate delayed terms, corresponding segment indicies, 
%and interpolation coefficients in the multi point boundary value
%problem of periodic orbits
% Input:
%   u0: state vector (N*M*n)
%   Ti: segment lengths(N)
%   M: dimension of used Chebyshev mesh
%   tau: time delay array
%   res: output resolution on each segment (default M)
% Output:
%   ud: state at t-tau(k) (M*n*N x n x nt)
%   id: segment index of t-tau(k) mapped back between 0 and T (M*n*N x nt)
%   fi: lagrange interpolation coefficients (M*n*N x nt x M)
%   nd: index of the periodic orbit where the tq is found without looping
%   around (0 current, 1,2,... past orbits) (M*n*N x nt)
%   dT: derivative of the querry point wrt Ti without looping around 
%   (M*n*N x nt x N x n_tau)

% Initialization
nt = length(tau);   % number of delays
N = length(Ti);     % number of segments
n = length(u0)/(N*M); % state dimension
ntau = ceil(max(tau)/sum(Ti)); % nt*T is neccessary to cover tau_max
if nargin < 5
    [ts,~] = bvp2sig(u0,Ti,M); % time mesh
else
    [ts,~] = bvp2sig(u0,Ti,M,res); % time mesh
end
ud = zeros(n,nt,length(ts)); % delayed terms (state at t_tau)
id = zeros(length(ts),nt); % index of the segment where delayed terms were interpolated
fi = zeros(length(ts),nt,M); % lagrange interpolation coefficients
nd = zeros(length(ts),nt); % index of present/past orbit where tq is found
dT = zeros(length(ts),nt,N,ntau+1); % derivative of tq wrt segment lengths

% Mapping points back to the periodic orbit
for i = 1:length(ts)
    
    % Span all possible delays
    for j = 1:length(tau)
        
        % Current segment index
        if nargin<5
            i0 = floor((i-1)/M)+1;
        else
            i0 = floor((i-1)/res)+1;
        end
        
        % Initialize derivatives
        delta_Ti = zeros(N,ntau+1); % wrt Ti
    
        % Coefficients of Ti in delta_i
        delta_i = ts(i)-sum(Ti(1:i0-1)); % distance from the starting point
        delta_Ti(i0,end) = delta_i/Ti(i0);  % account for initial segment
    
        % Step back segment by segment until tau is covered
        while delta_i<tau(j)
            % Cirle around segments
            if i0>1
                i0 = i0-1;
            else
                i0 = N;
                nd(i,j) = nd(i,j)+1;
            end
            
            % Emergency switch
            if nd(i,j)>ntau
                warning('Maximum iteration number reached while searchig for t_tau');
                break
            end
            
            % Update delta and its derivatives
            delta_i = delta_i + Ti(i0);
            delta_Ti(i0,end-nd(i,j)) = delta_Ti(i0,end-nd(i,j)) + 1;
            
        end
        
        % Identify interpolation segment and querry point
        id(i,j) = i0;
        tq = (delta_i-tau(j))/Ti(i0);
        ui0 = u0((i0-1)*M*n+1:i0*M*n);
        
        % Lagrange interpolation
        fiq = lag_coeff(M,tq);
        fi(i,j,:) = fiq;
        ud(:,j,i) = (fiq*reshape(ui0,M,n)).';
        
        % Evaluate derivative of t_tau wrt Ti
        dT(i,j,:,:) = 1/Ti(i0)*delta_Ti;
        dT(i,j,i0,end-nd(i,j)) = -(delta_i-tau(j)-...
            delta_Ti(i0,end-nd(i,j))*Ti(i0))/Ti(i0)^2;
        
    end
    
end

end

