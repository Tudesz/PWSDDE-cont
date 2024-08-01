function [x1,err,it] = newton_iter(x0,func,opts)
%NEWTON_ITER Newton iteration routine for solving a system of nonlinear
%equations
% Input:
%   x0: initial variable vector
%   func: a function, which returns the NAE and its Jacobian: [F, JF]
%      solving for F(x0) = 0, with JF(x0) = df(x)/dx at x0
%   opts: solver options structure
%    -> logs: if true print progress at each iteration step (default false)
%    -> reltol: stopping condition on correction step norm (default 1e-7)
%    -> abstol: stopping condition on error norm (default 1e-10)
%    -> maxiter: maximum number of iteration steps (default 10)
%    -> plots: if true plot progress of error function (default false)
% Output:
%   x1: solution of nonlinear equation
%   err: error of last funtion evaluation
%   it: number of iteration steps used

% Solver default options
if nargin < 3
    opts = corr_opts();
    opts = opts.nr;
end

% Initialize logs and plots
% if opts.logs
%     fprintf('\nNewton iteration\n')
% end
if opts.plots
    figure('Name','Newton iteration'); 
    subplot(2,1,1);
    xlabel('$i$'); ylabel('$|F_i|$'); 
    hold on; box on;
    set(gca, 'YScale', 'log')
    subplot(2,1,2);
    xlabel('$i$'); ylabel('$\Delta x_i$'); 
    hold on; box on;
end

% Run iteration
x1 = x0;
for it = 1:opts.maxiter
    tic;
    % Solution correction
    [err,JF] = func(x1);
    dx = -JF\err;
    x1 = x1 + dx;
    
    % Plot and log progress
    if opts.logs
        fprintf('   Iter %i, |dx|=%d, |F(x)|=%d, time %0.3f s\n',...
            it,norm(dx),norm(err),toc);
    end
    if opts.plots
        color = ones(1,3)*(opts.maxiter-it)/opts.maxiter;
        subplot(2,1,1)
        plot(abs(err),'.-','Color',color);
        title(sprintf('$||F(x_{%i})||$ = %d',it,norm(err)));
        subplot(2,1,2)
        plot(dx,'.-','Color',color);
        title(sprintf('$||\\Delta x_{%i}||$ = %d',it,norm(dx)));
        drawnow
        pause(1.0)
    end
    
    % Stop condition
    if norm(dx)<opts.reltol && norm(err)<opts.abstol
        break
    end
    if it==opts.maxiter && norm(err)>1e3*opts.abstol
        warning('Maximum iteration number reached with error %d',norm(err)); 
    end
end
if opts.plots
    hold off;
end

end

