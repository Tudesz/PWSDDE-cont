function plot_br1_ampl(branch,p_i,u_i,annotate,ls)
%PLOT_BR1_AMPL Display continuation result amplitudes in a state 
% dimension wrt 1 continuation parameter
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> p:  bifurcation parameters
%    -> sig: solution signature (n)
%   p_i: bifurcation parameter index
%   u_i: index of the state dimension whose amplitudes are displayed
%   annotate: level of annotations on the bifurcation plot (default 1)
%     -> 0: no extra information (ls options are applied)
%     -> 1: display stability and bifurcation points
%     -> 2: display only bifurcation points
%     -> 3: display only stability
%   ls: linestyle data forwarded to the plot function (default empty)

if nargin<4
    annotate = 1; % display stability and bifurcation points
end
if nargin<5 
    ls = {}; % no listyle data
elseif ~iscell(ls)
    ls = {ls};
end

% Preprocess output data
[p,ampl,~,mu_c,ibif] = get_br_data(branch,p_i(1),u_i);

% Plot results
switch annotate
    case 0
        plot(p,ampl,ls{:});
    case 1
        pi = p; pi(abs(mu_c)<=1) = NaN;
        plot(p,ampl,'b',pi,ampl,'r'); hold on
        mark_br_bifs(ibif,[p; ampl]);
    case 2
        plot(p,ampl,'k'); hold on;
        mark_br_bifs(ibif,[p; ampl]);
    case 3
        pi = p; pi(abs(mu_c)<=1) = NaN;
        plot(p,ampl,'b',pi,ampl,'r');
end
xlabel(sprintf('$\\lambda_%i$',p_i)); 
ylabel(sprintf('$(\\max(u_%i) - \\min(u_%i))/2$',u_i,u_i));
box on;

end

