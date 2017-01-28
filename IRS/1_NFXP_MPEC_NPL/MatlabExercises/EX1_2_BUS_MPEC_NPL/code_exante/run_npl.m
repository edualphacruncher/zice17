close all
clear all
clc
global ev0;  ev0 = 0;

%% Setup options and parameters.
% NFXP options
nfxp_options = nfxp.setoptions; % Set default options
nfxp_options.printfxp = 0;      % Print iteration info for fixed point algorithm if > 0

% Swiches:
bustypes = 4; % Select bus types smaller than this number: choose 3 or 4

% Spaces
mp.n    = 175;         % Number of gridpoints        
mp.max  = 450;         % Max of mileage
mp.grid = (0:mp.n-1)'; % Grid over mileage
mp.beta = 0.9999;      % Discount factor
mp.title = 'Nested Pseudo Likelihood (NPL)';

% Structural parameters for the policy iterator and FXP solution methods
% This is just an example to show that Psi (fixed-point equation in CCPspace)
%and Gamma (usual fixed-point equation in V space) produce  the same solution
% even though the methods are quite different, i.e. we don't estimate any parameters, just apply the
% Psi and Gamma mappings to update our CCPs and get the results.
mp_solve = mp;
mp_solve.RC = 9.7686;     % Replacement cost
mp_solve.c  = 1.3428;     % Cost parameter
mp_solve.p  = [0.0937,... % pi_0
               0.4475,... % pi_1
               0.4459,... % pi_2
               0.0127]';  % pi_3

% Create transition matrix
P_solve = nfxp.statetransition(mp_solve.p, mp_solve.n); %non-parametrically

% Solve the model using the policy iterator mapping. 
pk_npl  = npl.solve(mp_solve, P_solve);


% Solve the model using the hybrid FXP algorithm. pk_fxp are the CCPs using fxp.
% this is called FXP and not NFXP because we don't use the outer MLE loop (we assume (and are right)
% that our parameters have converged, so now we just need to find the solution to our model (CCPs)
% for these parameters.
[~,pk_fxp] = nfxp.solve(ev0,P_solve,0.001*mp_solve.c*mp_solve.grid,mp_solve,nfxp_options);

figure(1)
hold all
plot(mp.grid, 1-pk_npl,'LineWidth', 2);
plot(mp.grid,1-pk_fxp,'r', 'LineWidth', 1.5);
legend('\Psi', 'FXP', 'Location', 'southeast')
keyboard

% ************************************
% 0: Load data 
% ************************************

data = nfxp.readbusdata(mp, bustypes);

% ************************************
% 1: ESTIMATE p 
% ************************************
% Step 1) Set K=0
% Step 2) Estimate theta_p
tab = tabulate(data.dx1); % If you do not know what this does, print it
tab = tab(tab(:,3)>0,:);  % As above
m.p = tab(1:end-1,3)/100; % As above
P   = nfxp.statetransition(m.p, mp.n); % Create the transition matrix based on pi_0, pi_1, pi_2, and pi_3. Sparse.

% Step 3) Guess theta_c and CCPs (pk0)
mp.RC = 0; % Replacement cost
mp.c  = 0; % Maintenance parameter
theta_start = [mp.RC;mp.c]; % Make theta's for first likelihood estimation
pk0 = ones(mp.n,1)*0.01;    % Starting values for CCPs

% Step 4)  Maximize the likelihood
%       a) Construct the variables needed to calculate V_sigma
v = npl.phi(mp, pk0, P); % Calculate the variables used: pv_z(K), pv_z(R), pv_e 

%       b) Maximize the likelihood function given CCPs (still K=0)
options =  optimset('Display','off', 'GradObj','on', 'TolFun',1E-8,'TolX',0,'Hessian','on');
[theta_npl,FVAL,EXITFLAG,OUTPUT1] = fminunc(@(theta) npl.clogit(data, v, mp, theta),theta_start,options);

% Step 5) Update CCPs, remember to update given old CCPs, and thus old "v" (CCPs come from Logit=Psi(phi(pk)))
% phi)
pk = npl.pk(mp,theta_npl, v);

% Ignore step 6; if we stopped now, this would be Hotz-Miller!

for K = 1:100
    % Step 4a) again for updated pk. This is the phi mapping step.
    v  = npl.phi(mp, pk, P);    
    theta_old = theta_npl;

    % Step 4b) again for updated v (=phi(pk)) to update theta_npl using the old theta_npl as starting
    % value. This is the Psi mapping step
    [theta_npl,FVAL,EXITFLAG,OUTPUT1] = fminunc(@(theta) npl.clogit(data, v, mp, theta),theta_npl,options);
    
    % Step 5) again using updated theta_npl and old v (pk) 
    pk = npl.pk(mp,theta_npl, v);
    mp.RC = theta_npl(1);
    mp.c  = theta_npl(2);

    % Step 6)
    NPL_metric = max(abs(theta_old-theta_npl));
    mec_vec(K) = NPL_metric;

    mp.llval = npl.clogit(data, v, mp, theta_npl); % just get the likelihood value
    fprintf('K = %-10d %14.3f %14.3f %14.1f \n ' , ...
    K, ...
    mean(mp.RC) , ...
    mean(mp.c) , ...
    mean(mp.llval)*size(data.x,1)); 

    if NPL_metric < 1e-6
        break
    end
end

% Save likelihood value
mp.llval = npl.clogit(data, v, mp, theta_start);

%% Plot results
% Plot CCPs
figure(2)
plot(mp.grid,1-pk,'-r');
title('Choice probabilities ')
xlabel('Milage grid')
ylabel('Replacement probability')

% Plot convergence metric 
figure(3)
plot(2:K,mec_vec(2:end));
title('Convergence measure')
ylabel('supnorm(\theta_{k+1}+\theta_k)')
xlabel('Major iteration')


%% Print estimation results
fprintf('*************************************************************************\n'); 
fprintf('Fixed parameters and spaces\n'); 
fprintf('*************************************************************************\n'); 
fprintf('Beta           = %10.5f \n',mp.beta);
fprintf('n              = %10.5f \n',mp.n);
fprintf('Sample size    = %10.5f \n',size(data.x,1));
fprintf('\n'); 
fprintf('*************************************************************************\n'); 
fprintf('Method, %s\n', mp.title); 
fprintf('*************************************************************************\n'); 
fprintf('\n'); 
fprintf('Parameter Estimates\n');
fprintf('---------------------------------------------------------------------\n');
fprintf('%14s %14s  %14s %14s\n', '', 'RC', 'C' ,'Likelihood');
fprintf('---------------------------------------------------------------------\n');
fprintf('%14s %14.3f %14.3f %14.1f \n', 'param', ...
mean(mp.RC) , ...
mean(mp.c) , ...
mean(mp.llval)); 
fprintf('---------------------------------------------------------------------\n');
fprintf('\n\n');
