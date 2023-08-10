function [px,py] = p_bldosc(t,x,p)
%P_BLDOSC Periodic orbit plot functions for delayed bilinear oscillator
% Input:
%   t: time mesh of periodic orbit
%   x: state vector on t
%   p: system parameter vector
% Output:
%   px: x coordinates to be plotted
%   py: y coordinates to be plotted

px = x(1,:);
py = x(2,:);

end