function anim_br_orb(branch,sys,fs,ls,filename)
%ANIM_BR_ORB Animate orbit shape progression along a continuation output 
% branch
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> sig: solution signature (n)
%    -> p: paramter vector
%   sys: names of the functions that define the system
%    -> tau: time delay and its parameter Jacobians
%    -> tau_no: number of distinct time delays
%   fs: functions for plotting, default: fs1(t,y)=t wrt fs2(t,y)=y
%   ls: linestyle data forwarded to the plot function (default empty)
%   filename: if given save .avi file with this name
% Output:
%   video file with the name filename.avi

if nargin<4 || isempty(ls)
    ls = {}; % no listyle data
end

% Save as video file
if nargin>4
    myVideo = VideoWriter(filename,'Motion JPEG AVI');
    myVideo.Quality = 95; 
    myVideo.FrameRate = 10;
    open(myVideo)
end

figure();
for i = 1:length(branch)
    % Plot each orbit one by one
    if nargin>2 && ~isempty(fs)
        plot_orb(branch(i),sys,fs,100,ls)
    else
        plot_orb(branch(i),sys,[],100,ls)
    end
    title(sprintf('Orbit: %i',i));
    drawnow
    % Save frame to video file
    if nargin>3
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    else
        pause(0.1);
    end
end

if nargin>3
    close(myVideo);
end

end

