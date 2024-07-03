function anim_br_spectrum(branch,filename)
%ANIM_BR_ORB Animate orbit spectrum progression along a continuation output 
% branch
% Input:
%   branch: continuation run output structure
%    -> u: solution vectors (M*n)
%    -> T: segment lengths (N)
%    -> sig: solution signature (n)
%    -> p: paramter vector
%   filename: if given save .avi file with this name
% Output:
%   video file with the name filename.avi

% Save as video file
if nargin>1
    myVideo = VideoWriter(filename,'Motion JPEG AVI');
    myVideo.Quality = 95; 
    myVideo.FrameRate = 10;
    open(myVideo)
end

figure();
for i = 1:length(branch)
    % Plot the spectrum of each orbit one by one
    plot_spectrum(branch(i).mu)
    title(sprintf('Orbit: %i',i));
    drawnow
    % Save frame to video file
    if nargin>1
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    else
        pause(0.1);
    end
end

if nargin>1
    close(myVideo);
end

end

