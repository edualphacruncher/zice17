%**************************************************************************
% Projection Methods illustrated by solving simple Ramsey model
%**************************************************************************
% This script solves the deterministic Ramsey model using three different
% projection methods: 
% 1) Piecewise linear collocation
% 2) Chebyshev collocation
% 3) Chebyshev Galerkin (left as excercise)
%**************************************************************************
% Uses PlCollEqu.m and ChebCollEqu.m
%**************************************************************************
% By Johannes Brumm, 02/2016
%**************************************************************************

%% Initialization:

clear all;
tic % start clock

% Economic paramters to be chosen: 
beta  = 0.95; % discount factor
gamma = 2;    % risk aversion / inverse of intertemporal elasticity
delta = 0.1;  % rate of depreciation
alpha = 0.3;  % capital share

% Computational parameters to be chosen:
nrgps  = 20; % number of gridpoints
method_choice =1; % choose method: 
% 1=PL Collocation, 2=Chebyshev Collocation, 3=Chebyshev Galerkin
plot_choice = 0; % 0=make no plots; 1=make plots
op=optimset('Display','off','TolFun',10^-10,'MaxIter',100,'TolX',10^-10);
% set options for solver

% Parameters to be determined:
kstar=((1/beta-(1-delta))/alpha)^(1/(alpha-1)); % steady state capital stock
k_min = kstar*0.1; % minimal capital stock
k_max = kstar*1.5; % maximal capital stock
if method_choice==1 % case of piecewise linear interpolation
    gridstepsize = (k_max-k_min)/(nrgps-1);  % distance between two consecutive grid points
    kgrid = linspace(k_min,k_max,nrgps); % create grid between k_min and k_max with nrgps gridpoints
    kpgrid_start = 0.05*mean(kgrid)*ones(1,nrgps)+0.95*kgrid; % make starting guess for policy function
else % case of Chebyshev interpolatin 
    cheb_zeros=-cos((2*(1:nrgps)-1)*pi/(2*nrgps)); % calculate Chebyshev zeros
    kgrid=(cheb_zeros+1)*(k_max-k_min)/2+k_min; % create Chebyshev grid
    cheb_M=cos((0:1:nrgps-1)'*acos(2*(kgrid-k_min)/(k_max-k_min)-1)); % Matrix Chebyshev polynomials (row) evaluated at Chebyshev zeros (column)
    kpgrid_start = 0.05*mean(kgrid)*ones(1,nrgps)+0.95*kgrid; % simple initial guess for the capital-policy at gridpoints
    cheb_co_start=(cheb_M*kpgrid_start')./sum((cheb_M').^2)'; % translate guess for capital-policy into guess for Chebyshev coefficients
end
Errorgrid=linspace(k_min,k_max,1000); % define a fine grid for error analysis


%% Main step (Solve for parameters of policy function):

if method_choice==1
    method='Piecewise Linear Collocation';
    [kpgrid,f,EXITFLAG,OUTPUT,JACOB]=fsolve(@(v) PlCollEqu(v,kgrid,nrgps,gridstepsize,beta,gamma,delta,alpha),kpgrid_start,op);
elseif method_choice==2
    method='Chebyshev Collocation';
    [cheb_co,f,EXITFLAG,OUTPUT,JACOB]=fsolve(@(v) ChebCollEqu(v,kgrid,nrgps,k_min,k_max,beta,gamma,delta,alpha),cheb_co_start,op);
elseif method_choice==3
    method='Chebyshev Galerkin';
    % This is left as an excercise!
end
toc % stop clock

%% Error Analysis:

EulerErrors=zeros(1,size(Errorgrid,2)); % initialize Euler error vector
Policygrid=Errorgrid; % initialize a grid for the policy (needed for plot) 
i=0;
for k = Errorgrid % go through the entire error grid
    i=i+1;
    
    if method_choice==1
    kp = interp1(kgrid,kpgrid,k); % determine capital choice this period
    kpp = interp1(kgrid,kpgrid,kp); % determine capital choice next period
    else
    kp = cheb_co'*cos((0:nrgps-1)'*acos(2*(k-k_min)/(k_max-k_min)-1)); 
    % determine capital choice this period by Chebyshev interpolation
    kpp = cheb_co'*cos((0:nrgps-1)'*acos(2*(kp-k_min)/(k_max-k_min)-1)); 
    % determine capital choice next period by Chebyshev interpolation
    end
    Policygrid(i)=kp; % save policy 
    EulerErrors(i) = abs(((k^alpha + k*(1-delta) - kp) / (beta*(alpha*kp^(alpha-1) + (1-delta))*(kp^alpha + kp*(1-delta)- kpp)^(-gamma))^(-1/gamma)) - 1);
    % calculate Euler error at the point on the error grid
end

MaxEulerError=max(EulerErrors); % calculate the maximum Euler error 
AverageEulerError=mean(EulerErrors); % calculate the average Euler error 
display(['Method: ',method,'. ',num2str(nrgps),' Points.'])  
display(['Maximum EE: ',num2str(MaxEulerError),'. Average EE: ',num2str(AverageEulerError),'.'])


%% Plot Results:

if plot_choice==1

% Plot Euler errors:
figure
plot(Errorgrid,EulerErrors)
legend(['Time-Iteration with ', num2str(nrgps),' Grid Points'],'Location','northwest')
title(['~~~~Det.~Ramsey: $\gamma=$', num2str(gamma),', $\beta=$',num2str(beta),', $\alpha=$',num2str(alpha),', $\delta=$',num2str(delta)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Capital stock $k$','interpreter','latex','FontSize',12)
ylabel('Euler error $E(k)$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')

% Plot policy function:
figure
plot(Errorgrid,Policygrid)
legend(['Time-Iteration with ', num2str(nrgps),' Grid Points'],'Location','northwest')
title(['Det.~Ramsey: $\gamma=$', num2str(gamma),', $\beta=$',num2str(beta),', $\alpha=$',num2str(alpha),', $\delta=$',num2str(delta)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Capital stock $k_{t}$','interpreter','latex','FontSize',12)
ylabel('Capital choice $k_{t+1}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
text(kstar,kstar,' \leftarrow steady state')

end

