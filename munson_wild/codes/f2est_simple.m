function [fder2,s2n] = f2est_simple(evalf,nf,xb,h,p,fval,fnoise,varargin)
%
%   [fder2,s2n] = f2est_simple(evalf,nf,xb,h,p,fval,fnoise,varargin)
%
% This is an experimental code segment for estimating the second
% derivative. The f7pt value may be a better estimate than
% the mean of the difference provided all the 2nd order
% differences have the same sign.
%
%     Jorge More' and Stefan Wild. Argonne National Laboratory, 2010.
%
% Based on the paper
%   Estimating Derivatives of Noisy Simulations, ACM Transactions on 
%   Mathematical Software, 38(3):19:1-19:21, 2012. DOI: 10.1145/2168773.2168777
% See http://www.mcs.anl.gov/~wild/cnoise for additional details/references

mid = (length(fval) + 1)/2;
fc = fval(mid);
  
  % Estimate the second derivative with fnoise^(1/4).

  hval = fnoise^(0.25);
  xp = xb + hval*p;
  fr = feval(evalf,xp,varargin{:});
  xp = xb - hval*p;
  fl = feval(evalf,xp,varargin{:});
  f3pt_a = ((fr-2*fc+fl)/hval)/hval; % Recall that fc was saved earlier.

  if (fnoise > 0); s2n = abs(fr-2*fc+fl)/fnoise; else s2n = 0; end
  fmax = max(abs([fr,fc,fl]));
  fmin = min(abs([fr,fc,fl]));
  ftest = (fmax - fmin)/fmax;
  
  fder2 = 0;
  if (s2n >= 100 && ftest <= 0.1) 
    fder2 = f3pt_a;
  else
    ht = (fnoise/abs(f3pt_a))^(1/4);
%    if (s2n >= 100 && ftest <= 0.5 && (ht >= hval/10 && ht <= 10*hval))
%      fder2 = f3pt;
%    end
  end

%   fprintf('\n\n%s %9.2e %9.2e %9.2e %9.2e\n',' h, fder2, s2n ',...
%           hval,f3pt_a,s2n,ftest);

  if (fder2 == 0)
    hval = ht;
    xp = xb + hval*p;
    fr = feval(evalf,xp,varargin{:});
    xp = xb - hval*p;
    fl = feval(evalf,xp,varargin{:});
    f3pt_b = ((fr-2*fc+fl)/hval)/hval; % Recall that fc was saved earlier.
    
    if (fnoise > 0); s2n = abs(fr-2*fc+fl)/fnoise; else s2n = 0; end
    fmax = max(abs([fr,fc,fl]));
    fmin = min(abs([fr,fc,fl]));
    ftest = (fmax - fmin)/fmax;
    if (s2n >= 100 && ftest <= 0.1) 
      fder2 = f3pt_b;
    end
    if (abs(f3pt_a-f3pt_b) <= 0.5*max(abs([f3pt_a,f3pt_b])))
      fder2 = f3pt_b;
    end
%     fprintf('%s %9.2e %9.2e %9.2e %9.2e\n\n',' h, fder2, s2n ',...
%             hval,f3pt_b,s2n,ftest);
  end

%   if (fder2 == 0)
%       ht = (fnoise/abs(f3pt_b))^(1/4);
%     hval = ht;
%     xp = xb + hval*p;
%     fr = feval(evalf,xp,varargin{:});
%     xp = xb - hval*p;
%     fl = feval(evalf,xp,varargin{:});
%     f3pt_c = ((fr-2*fc+fl)/hval)/hval; % Recall that fc was saved earlier.
%     
%     if (fnoise > 0); s2n = abs(fr-2*fc+fl)/fnoise; else s2n = 0; end
%     fmax = max(abs([fr,fc,fl]));
%     fmin = min(abs([fr,fc,fl]));
%     ftest = (fmax - fmin)/fmax;
%     if (s2n >= 100 && ftest <= 0.1) 
%       fder2 = f3pt_c;
%     end
%     if (abs(f3pt_c-f3pt_b) <= 0.5*max(abs([f3pt_c,f3pt_b])))
%       fder2 = f3pt_c;
%     end
%     fprintf('%s %9.2e %9.2e %9.2e %9.2e\n\n',' h, fder2, s2n ',...
%             hval,f3pt_b,s2n,ftest);
%   end
% 
%   if fder2==0
% fder2 = median([ f3pt_a f3pt_b f3pt_c])
%      pause
%   end