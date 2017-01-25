% This script runs tests on a simple n-dimensional noisy quadratic. 
% Stefan Wild. Argonne National Laboratory, January 2016.

rand('state',1); % Matlab may warn, but here's how I get reproducibility

% Define the function we're going to consider
func=@examplefun; 
global probnum L sigma  % This is used to define functions in examplefun
probnum = 1; 
n = 3;          % Number of dimensions
x = rand(n,1);  % The point of interest
p = rand(n,1);  % The direction along which to compute derivative
fderact = 2*x'*p; % "True" directional derivative (make empty for your function)
if probnum==3
    fderact = 2*( x'.*10.^(2*(1-(1:n))) )*p;
end
 
 
% Step 1: 
% Evaluate the function slong the direction p
h = 1e-10;
nf = 9;
for i=1:nf
    fval(i) = func(x + h*(i-5)*p);
end

figure(1)
plot(h*[-4:4],fval,'-s')
ylabel('f(x + u*p)'), xlabel('u'), title('Function values')


% Step 2: 
% Try a whole bunch of forward-difference parameter values to f'
hvec=logspace(-12,0,40); fderh=zeros(size(hvec));
for j=1:length(hvec)
    fderh(j) = (func(x + hvec(j)*p)-fval(5))/hvec(j);
end
figure(2) % !! Note: this is only justified for f'(x)>0
loglog(hvec,abs(fderh-fderact)/abs(fderact),'x');
ylabel('| f''(x) - ( f(x+h*p)-f(x) ) / h | / |f''(x)|'), xlabel('h'), 
title('Relative error of forward-difference estimate')


% Step 3:
% Compute noise estimate
[fnoise,level,inform] = ECnoise(nf,fval);

if (inform == 1) infos = 'Noise detected';
    fprintf([' ',infos,'      %9.2e\n'],fnoise);
elseif (inform == 2) infos = '*Noise not detected; h too small'; disp(infos); return
elseif (inform == 3) infos = '*Noise not detected; h too large'; disp(infos); return
end


% Step 4:
% Get a coarse estimate of f"  (alternatively, set fder2=max(1,abs(f(x)));) 
% Estimate derivative
[fder2,s2n] = f2est_simple(func,nf,x,h,p,fval,fnoise);
if (s2n >= 50) 
    fprintf(' Estimated fder2     %9.2e  (s2n=%9.2e)\n',fder2, s2n)
else
    fprintf('*Procced with caution, low signal to noisein fder2!*\n')
    fprintf(' Estimated fder2     %9.2e  (s2n=%9.2e)\n',fder2, s2n)
end    

hopt = 1.68*sqrt(fnoise/abs(fder2));
fderest = (func(x + hopt*p)-fval(5))/hopt;
fprintf(' Forward-diff step   %9.2e\n',hopt)
fprintf(' Forward diff est    %9.2e  ("true"=%9.2e)\n',fderest, fderact)

figure(2)
hold on 
plot(hopt,abs(fderest-fderact)/abs(fderact),'go');






    
    