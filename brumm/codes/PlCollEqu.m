function [error_EE]=PlCollEqu(v,kgrid,nrgps,gridstepsize,beta,gamma,delta,alpha)

kpgrid=v;

   for gp = 1:nrgps % go through entire grid       
       k = kgrid(gp); % read grid point to be considered
       kp = kpgrid(gp); % read capital choice to be considered
       % determine kpp by linear spline interpolation:
       gp_right = find(kgrid>kp,1);  
       kpp = ((kgrid(gp_right)-kp)*kpgrid(gp_right-1) + (kp-kgrid(gp_right-1))*kpgrid(gp_right))/gridstepsize;
       % alternative but slower interpolation: kpp = interp1(kgrid,kpgrid,kp);
       error_EE(gp)=(k^alpha + k*(1-delta) - kp)^(-gamma) - beta*(alpha*kp^(alpha-1) + (1-delta))*(kp^alpha + kp*(1-delta)- kpp)^(-gamma);
   end  
end

