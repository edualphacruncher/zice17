%**************************************************************************
% Time iteration illustrated by solving simple Ramsey model
%**************************************************************************
% This script solves the deterministic Ramsey model using a simple 
% time iteration algorithm. 
% The Euler equation is solved at each grid point using Newton's method.
% The policy function is interpolated by a linear spline.
%**************************************************************************
% By Johannes Brumm, 01/2016
%**************************************************************************


%% Initialization:
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
gridstepsize = (k_max-k_min)/(nrgps-1); % distance between two consecutive grid points
kgrid = linspace(k_min,k_max,nrgps); % create grid between k_min and k_max with nrgps gridpoints
kpgrid = 0.05*mean(kgrid)*ones(1,nrgps)+0.95*kgrid; % make starting guess for policy function
error_it = 1; % initial value for error in timeiteration
Errorgrid=linspace(k_min,k_max,1000); % define a fine grid for error analysis


%% Time Iteration:

% Run time-iteration algorithm until error below tolerance
while error_it > tol_it
    
    for gp = 1:nrgps % loop over all gridpoints
       
       k = kgrid(gp); % read grid point to be considered
       kp = k; % use the current capital level as starting guess
       
       % Run Newton's Method to find solution at given gridpoint.
       error_EE=1;
       while abs(error_EE) > tol_EE 
           
           % Determine capital choice next period by linear interpolation:
             kpp = interp1(kgrid,kpgrid,kp); % given capital choice today (kp) and policy function from last time-iteration step (kpgrid) 
            
           % Check whether consumption today and tomorrow is positive:
           if (k^alpha + k*(1-delta) - kp) > 0 && (kp^alpha + kp*(1-delta)- kpp) > 0 
           
               % Calculate residual of Euler equation:
               error_EE=(k^alpha + k*(1-delta) - kp)^(-gamma) - beta*(alpha*kp^(alpha-1) + (1-delta))*(kp^alpha + kp*(1-delta)- kpp)^(-gamma);
               
               % Calculate (right) derivative of residual of Euler equation at kp:
               kpdelta=kp+10^-5;
               kppdelta = interp1(kgrid,kpgrid,kpdelta);
               error_EEdelta=(k^alpha + k*(1-delta) - kpdelta)^(-gamma) - beta*(alpha*kpdelta^(alpha-1) + (1-delta))*(kpdelta^alpha + kpdelta*(1-delta)- kppdelta)^(-gamma);
               derivative=(error_EEdelta-error_EE)/(10^-5);
               kp=min(k_max-2*10^-5,max(k_min,kp-(derivative^(-1))*error_EE));
               
           elseif (k^alpha + k*(1-delta) - kp) < 0
               kp = kp*0.99; % move guess to the left, if consumption today is below 0 
           elseif (kp^alpha + kp*(1-delta)- kpp) < 0
               kp = kp*1.01; % move guess to the right, if consumption tomorrow is below 0          
           end
       end
       kpgridnew(gp)=kp; % save new policy    
    end
    % Calculate error:
    error_it=max(abs(kpgrid-kpgridnew))/(max(abs(kpgrid)));
    % Update policy:
    kpgrid=kpgridnew;

end
toc % stop clock


%% Error Analysis:

EulerErrors=zeros(1,size(Errorgrid,2)); %initialize Euler error vector
i=0;
for k = Errorgrid % go through the entire error grid
    i=i+1;
    kp = interp1(kgrid,kpgrid,k); % determine capital choice today by linear interpolation
    kpp = interp1(kgrid,kpgrid,kp); % determine capital choice next period by linear interpolation
    % Calculate Euler equation error:
    EulerErrors(i) = abs(((k^alpha + k*(1-delta) - kp) / (beta*(alpha*kp^(alpha-1) + (1-delta))*(kp^alpha + kp*(1-delta)- kpp)^(-gamma))^(-1/gamma)) - 1);    
end

MaxEulerError=max(EulerErrors); % calculate the maximum Euler error 
AverageEulerError=mean(EulerErrors); % calculate the average Euler error 
display(['Method: Time Iteration With Linear Splines. ',num2str(nrgps),' Points.']) 
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
plot(kgrid,kpgrid)
legend(['Time-Iteration with ', num2str(nrgps),' Grid Points'],'Location','northwest')
title(['Det.~Ramsey: $\gamma=$', num2str(gamma),', $\beta=$',num2str(beta),', $\alpha=$',num2str(alpha),', $\delta=$',num2str(delta)],'interpreter','latex')
set(gca,'FontSize',12)
xlabel('Capital stock $k_{t}$','interpreter','latex','FontSize',12)
ylabel('Capital choice $k_{t+1}$','interpreter','latex','FontSize',12)
set(legend,'FontSize',12,'interpreter','latex')
text(kstar,kstar,' \leftarrow steady state')

end