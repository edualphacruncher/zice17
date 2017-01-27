clear;
global ev0;  ev0=0;

%% NFXP options
nfxp_options=nfxp.setoptions;  			    % Set default options
nfxp_options.printfxp   = 1            	% Print iteration info for fixed point algorithm if > 0

nfxp_options.min_cstp=100;										% Set miniumn number of contraction steps (successive approximations)
nfxp_options.max_cstp=nfxp_options.min_cstp;	% Set maximum number of contraction steps (successive approximations)

% Spaces
mp.n=175;				        % Number of gridpoints		
mp.max=450;				      % Max of mileage
mp.grid= (0:mp.n-1)';  	% Grid over mileage

% Structural parameters
mp.p=[0.0937 0.4475 0.4459 0.0127]';   	% Transition probabiliuties
mp.RC=11.7257;     						% Replacement cost
mp.c=2.45569;									% Cost parameter
mp.beta=0.9999;								% Discount factor

P0 = nfxp.statetransition(mp.p, mp.n);	% Transition matrix for mileage
cost0=0.001*mp.c*mp.grid;				% Cost function
tic
[mp.ev, pk0, F, iterinfo]=nfxp.solve(0, P0, cost0, mp, nfxp_options); 	% Solve the model
toc
