function y = examplefun(x)
%   function y = examplefun(x)
% This is an example of a simple n-dimensional noisy quadratic. 
% There are three optional global variables:
%   probnum [1] :   selects which version of the function is used
%                       1 = determinstic
%                       2 = stochastic
%                       3 = determinstic, unscaled
%  L [15] :         for probnum=1, increasing L generally increases noise
%  sigma [1e-4] :   for probnum=2, increasing sigma increases noise
% Stefan Wild. Argonne National Laboratory, January 2016.

global probnum L sigma
% Set default values:
if isempty(L), L = 15; end 
if isempty(probnum), probnum = 1; end
if isempty(sigma), sigma = 1e-4; end    
n = length(x);

    switch probnum
        case 1
            % This function is a fun way to compute x'*x
            y = norm(x)^2;
            for i=1:L, y = sqrt(y); end
            for i=1:L, y = y^2; end
        case 2
            % This is a stochastic function with
            %   E[f(x)] = x'*x
            % Var[f(x)] = (x'*x)^2 * sigma^2 
            y = norm(x)^2;
            y = y*(1+sigma*sqrt(3)*(2*rand-1));
        case 3
            % This function is a fun way to compute x'*x
            y = sum( (10.^(1-(1:n)) .* x(:)').^2 );
            for i=1:L, y = sqrt(y); end
            for i=1:L, y = y^2; end
    end
    