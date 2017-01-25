% This script runs the random gradient on a noisy quadratic. 
% It illustrates potential cost of using a zero-order method, role of
% difference parameter, effect of nonstationary noise, etc.
%
% Stefan Wild. Argonne National Laboratory, January 2016.


% Define problem and constants:
% (Don't forget probnum, L, sigma)
%-% try 
 global probnum L sigma
% probnum = 2; sigma=1e-9;
% probnum = 2; sigma=1e-18;
% func=@examplefun; 
 
n = 2;
maxit = 100;
delta = 1/(4*(n+4)*2);

for test=1:2 % Test different forward difference parameters
    if test==2
        h = 1e-11;
    elseif test==1
        h = 1e-5;
    end
    
    randn('state',1); % Matlab may warn, but here's how I get reproducibility
    
    % Starting point
    x = 10*randn(n,1);
    f(1) = func(x);
    
    for k=1:maxit
        d = randn(n,1); % Generate a Gaussian direction
        g = ( (func(x(:,k)+h*d) - f(k))/h ) * d; % Normalize step
        x(:,k+1) = x(:,k)-delta*g;
        f(k+1) = func(x(:,k+1));
        
        
        % Plotting
        figure(2+test); clf;
        subplot(1,2,1), semilogy(1:k+1,f(1:k+1)), xlim([1,maxit])
        ylabel('F value'), xlabel('Iterations'), legend(['k=',num2str(k)])
        subplot(1,2,2), plot(x(1,1:k+1),x(2,1:k+1),'x-',x(1,k+1),x(2,k+1),'go')
        ylabel('x_2 value'), xlabel('x_1 value')
        pause(0.1) % Just for vis
    end
    pause
end