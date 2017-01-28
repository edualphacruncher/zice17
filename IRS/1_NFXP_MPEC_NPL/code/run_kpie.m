% Strutural Estimation of Rust's engine replacement model using NPL 
close all
% clear all
clc

% Swiches:
Kmax=10 					% Max number of outer loop iterations for NPL 
bustypes = 4; % Select bus types smaller than this number: choose 3 or 4

% Spaces and parameters
mp.n    = 175;         % Number of gridpoints        
mp.max  = 450;         % Max of mileage
mp.grid = (0:mp.n-1)'; % Grid over mileage
mp.beta = 0.9999;      % Discount factor


% ************************************
% section 0: Load data 
% ************************************
data = nfxp.readbusdata(mp, bustypes);

% ************************************
% section 1: ESTIMATE p 
% ************************************
% Initialization: Estimate theta_p
tab = tabulate(data.dx1); % If you do not know what this does, print it
tab = tab(tab(:,3)>0,:);  % As above
m.p = tab(1:end-1,3)/100; % As above
P   = nfxp.statetransition(m.p, mp.n); % Create the transition matrix based on pi_0, pi_1, pi_2, and pi_3

% ************************************
% section : INITIAL CCP's 
% ************************************
% Initial guess for the conditional choice probabilities
if 0
	pk0 = ones(mp.n,1)*0.99;    % Starting values for CCPs
	% Staring values of structural parameters
	theta0 = 0*[9.7686;1.3428];  %  [mp.RC; mp.c]

% Ideally should be non-parametric estimates
else % Use static model to intialize NPL with CCPs
	mp0=mp;
	mp0.beta=0;
	theta0 = 0*[9.7686;1.3428];
	[mp0, pk0, logl0, K0]=npl.estim(theta0, pk0, data, P, mp0, 1);
	theta0=[mp0.RC, mp0.c];
end


% ************************************
% 2: Estimate Rust's model using NPL 
% ************************************
[mp, pk, logl, K]=npl.estim(theta0, pk0, data, P, mp, Kmax);

% *******************************************
% 3: Solve estimated model using NPL and FXP 
% *******************************************

% Solve for fixed point of policy iteration operator, pk=Psi(pk)
fignr =1; % Figure number for figer that displays convergence of Policy iteration operator. (optional argument in npl.solve)
pk_npl_fixpoint  = npl.solve(mp, P, pk0, fignr);  % Solve for  

% Solve for fixed point of bellman operator in EV space, pk=Psi(pk)the model using the hybrid FXP algorithm, 
nfxp_options = nfxp.setoptions; % Set default options

ev0=0;
[ev_fxp,pk_fxp] = nfxp.solve(ev0,P,0.001*mp.c*mp.grid,mp,nfxp_options);

% Plot CCPs
figure(2)
hold all
plot(mp.grid,1-pk,'-r','LineWidth', 2);
plot(mp.grid,1-pk_npl_fixpoint,'-b','LineWidth', 2);
plot(mp.grid,1-pk_fxp,'--k','LineWidth', 1.5);
title(sprintf('Replacement probability, K=%d', K));
legend('Last evaluation of \Psi','Fixed point of \Psi' ,'Fixed point of \Gamma','Location', 'southeast')
xlabel('Milage grid')
ylabel('Replacement probability')
grid on
ylim([0 0.16])

