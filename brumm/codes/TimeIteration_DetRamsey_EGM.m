%**************************************************************************
% Endogenous grid method (EGM) illustrated by solving simple Ramsey model
%**************************************************************************
% This script solves the deterministic Ramsey model 
% with wealth as state variable using a time iteration algorithm  
% with the endogenous grid method. 
% The policy function is interpolated by a linear spline.
%**************************************************************************
% By Johannes Brumm, 01/2016
%**************************************************************************

% Initialization:
clear all;
tic % start clock

% Economic paramters to be chosen: 
beta  = 0.95; % discount factor
gamma = 2; % risk aversion / inverse of intertemporal elasticity
delta = 0.1;  % rate of depreciation
alpha = 0.3;  % capital share

% Computational parameters to be chosen:
nrgps  = 20;    % number of grid points
tol_EE = 10^-8; % error tolerance for solution to the Euler equation
tol_it = 10^-6; % error tolerance of time-iteration algorithm
plot_choice = 0; % 0=make no plots; 1=make plots

% Parameters to be determined:
kstar=((1/beta-(1-delta))/alpha)^(1/(alpha-1)); % steady state capital stock
k_min = kstar*0.1; % minimal capital stock
k_max = kstar*1.5; % maximal capital stock
w_min = k_min^alpha + k_min*(1-delta); % minimal wealth
w_max = k_max^alpha + k_max*(1-delta); % maximal wealth
gridstepsize = (w_max-w_min)/(nrgps-1); % distance between two consecutive grid points
wgrid = linspace(w_min,w_max,nrgps); % create grud between w_min and w_max with nrgps gridpoints
%kpgrid = wgrid/2; % simple initial guess for the capital-policy function
kpgrid = linspace(k_min,k_max,nrgps); 
error_it = 1; % initial value for error in timeiteration
Errorgrid=linspace(w_min,w_max,1000); % define a fine grid for error analysis

%% Policy Function Iteration:

% Run time-iteration algorithm until error below tolerance
while error_it > tol_it
    
    for gp = 1:nrgps % loop over all gridpoints
       
        kp =  kpgrid(gp); % read grid point to be considered
           
        % Determine choice next period by linear interpolation:
        wp = kp^alpha+(1-delta)*kp;
        kpp = interp1(wgrid,kpgrid,wp);
           
        % get wealth w from choice kp
        w = kp + (beta*(alpha*kp^(alpha-1) + (1-delta))*(wp - kpp)^(-gamma))^(1/-gamma);
               
        wgridnew(gp)=w; % save new policy    
    end

    % Calculate error:
    error_it=max(abs(wgrid-wgridnew))/(max(abs(wgrid)));
    
    % Update policy:
    wgrid=wgridnew;
end
toc

%% Error Analysis:

EulerErrors=zeros(1,size(Errorgrid,2)); %initialize Euler error vector
i=0;
    for w = Errorgrid %go through the entire error grid
    i=i+1;
    kp = interp1(wgrid,kpgrid,w,'linear','extrap');
    wp = kp^alpha+(1-delta)*kp;
    kpp = interp1(wgrid,kpgrid,wp);
    EulerErrors(i) = abs(((w - kp) / (beta*(alpha*kp^(alpha-1) + (1-delta))*(wp- kpp)^(-gamma))^(-1/gamma)) - 1);   
end

MaxEulerError=max(EulerErrors); % calculate the maximum Euler error 
AverageEulerError=mean(EulerErrors); % calculate the average Euler error 
display(['Method: Time Iteration With EGM and Linear Splines. ',num2str(nrgps),' Points.']) 
display(['Maximum EE: ',num2str(MaxEulerError),'. Average EE: ',num2str(AverageEulerError),'.'])



%% Plot Results:

if plot_choice==1

% Plot Euler errors:
figure
plot(Errorgrid,EulerErrors)
legend(['Time-Iteration with EGM and ', num2str(nrgps),' Grid Points'],'Location','northwest')
title(['~~~~Det.~Ramsey: $\gamma=$', num2str(gamma),', $\beta=$',num2str(beta),', $\alpha=$',num2str(alpha),', $\delta=$',num2str(delta)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Wealth $w$','interpreter','latex','FontSize',12)
ylabel('Euler error $E(w)$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')

% Plot policy function:
figure
plot(wgrid,kpgrid)
legend(['Time-Iteration with EGM and ', num2str(nrgps),' Grid Points'],'Location','northwest')
title(['Det.~Ramsey: $\gamma=$', num2str(gamma),', $\beta=$',num2str(beta),', $\alpha=$',num2str(alpha),', $\delta=$',num2str(delta)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Wealth $w_{t}=f(k_t) + (1-\delta) k_t$','interpreter','latex','FontSize',12)
ylabel('Capital choice $k_{t+1}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')

end