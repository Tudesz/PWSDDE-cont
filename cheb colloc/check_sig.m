function sig_error = check_sig(orb,sys)
%CHECK_SIG Solution signature of non-smooth periodic orbit for continuity
% Input:
%   orb: starting periodic orbit data structure
%    -> sig: sloution signature (event list)
%   sys: names of the functions that define the 
%    -> e: event function, map and corresponding Jacobians
% Output:
%   sig_error: true if there is a continuity issue

% No error by default
sig_error = false;

% Check first event
if feval(sys.e,[],[],[],orb.sig(1),7,1) ~= ...
        feval(sys.e,[],[],[],orb.sig(end),7,2)
    sig_error = true;
end

% Loop through the rest of the events
for i = 2:length(orb.sig)
    if feval(sys.e,[],[],[],orb.sig(i-1),7,2) ~= ...
            feval(sys.e,[],[],[],orb.sig(i),7,1)
        sig_error = true;
    end
end

% Throw warning
if sig_error
    warning('Solution signature error!')
end

end

