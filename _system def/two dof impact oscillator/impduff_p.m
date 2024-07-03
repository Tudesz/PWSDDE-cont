function [px,py,xlab,ylab] = impduff_p(t,x,xd,p)
%IMPDUFF_P Periodic orbit plot functions for 2 DoF impact Duffing
%oscillator
% Input:
%   t: time mesh of periodic orbit
%   x: state vector on t
%   xd: delayed state vector on t
%   p: system parameter vector
% Output:
%   px: x coordinates to be plotted
%   py: y coordinates to be plotted
%   xlab: label on the x axis
%   ylab: label on the y axis


px = x(1,:);
py = x(2,:)- x(1,:);
xlab = '$x_1(t)$';
ylab = '$x_2(t)-x_1(t)$';
end

