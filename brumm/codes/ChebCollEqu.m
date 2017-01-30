function [error_EE]=ChebCollEqu(v,kgrid,nrgps,k_min,k_max,beta,gamma,delta,alpha)

cheb_co=v; 

% go through all gridpoints:
for gp = 1:nrgps        
    % read grid point to be considered and the respective capial choice:
    k = kgrid(gp); 
    % determine capital choice this period by Chebyshev interpolation:
    kp = cheb_co'*cos((0:nrgps-1)'*acos(2*(k-k_min)/(k_max-k_min)-1)); 
    % determine capital choice next period by Chebyshev interpolation:
    kpp = cheb_co'*cos((0:nrgps-1)'*acos(2*(kp-k_min)/(k_max-k_min)-1)); 
    % calculate consumption and ensure that it is positive:
    c=max(10^-6,(k^alpha + k*(1-delta) - kp)); cp=max(10^-6,(kp^alpha + kp*(1-delta)- kpp));
    % calculate error in Euler equations:
    error_EE(gp)=(c)^(-gamma) - beta*(alpha*kp^(alpha-1) + (1-delta))*(cp)^(-gamma);
end 

end

