global ev0;  ev0 = 0;

%% NFXP options
nfxp_options = nfxp.setoptions;  			% Set default options
nfxp_options.printfxp   = 0;            % Print iteration info for fixed point algorithm if > 0
% nfxp_options.hessian  = 'analytical';		

% Print nfxp options: 
nfxp_options = nfxp_options

% Swiches:
MonteCarlo = 1;   % 1: Do monte Carlo, 0: Use rust bus data 
twostep = 1;      % 1: Estimate using two-step pmle, 0: Use full mle (NOTE: MPEC NOT implemnted for twostep=0)
bustypes = 3;		% Select bus types smaller than this number, i.e. 1,2,..bustypes

% Spaces
mp.n = 175;             % Number of gridpoints		
mp.max = 450;          % Max of mileage
mp.grid = (0:mp.n-1)'; % Grid over mileage

% Parameters
mp.beta = 0.9999;                      % Discount factor

% Structural parameters initialized only for Monte Carlo
mp.p  = [0.0937 0.4475 0.4459 0.0127]';% Transition probabiliuties
mp.RC = 11.7257;                       % Replacement cost
mp.c  = 2.45569;                       % Cost parameter

if MonteCarlo   % Solve model for DGP
	nMC = 1; 	    % number of MC samples
	N = 50;	      % Number of busses to simulate 
	T = 119;		  % Number of time periods to simulate 

	P0 = nfxp.statetransition(mp.p, mp.n);	% Transition matrix for mileage
	cost0 = 0.001*mp.c*mp.grid;				% Cost function
	[mp.ev, pk0] = nfxp.solve(0, P0, cost0, mp, nfxp_options); 	% Solve the model
% ************************************************************************
fprintf('Begin Monte Carlo, with n = %d replications\n', nMC);
% ************************************************************************

else    % Read in and discretize Rust's data
	nMC = 1; 
 	timetosimulate = tic;
	data = nfxp.readbusdata(mp, bustypes);
	timetosimulate = toc(timetosimulate);
end

rand('seed',300);

%main loop
for i_mc = 1:nMC;
	% ************************************
	% STEP 0: SIMULATE DATA 
	% ************************************

	if MonteCarlo;
		timetosimulate = tic;
		data = nfxp.simdata(N, T, mp, P0, pk0);
		timetosimulate = toc(timetosimulate);
    fprintf('i_mc = %d, Time to simulate data : %1.5g seconds\n', i_mc, timetosimulate);
	end
		samplesize = numel(data.d);

	% ************************************
	% STEP 1: ESTIMATE p 
	% ************************************

	% Estimate model using Partial MLE
	tab = tabulate(data.dx1);
	tab = tab(tab(:,3)>0,:);
	p = tab(1:end-1,3)/100;
	P = nfxp.statetransition(p, mp.n);

  % ************************************
	% STEP 2a: ESTIMATE parameters using NFXP
	% ************************************

	% Starting values
	mp.beta = mp.beta;   
	mp.p = p;

  theta_start_pml = [0;0];
  ll_bhhh = @(theta)   nfxp.ll(data, mp, P, nfxp_options, theta);
    
	options =  optimset('Display','off', 'GradObj','on', 'TolFun',1E-8,'TolX',0,'Hessian','on');

  outsidetimer = tic;
	[theta_hat,FVAL,EXITFLAG,OUTPUT1] = fminunc(ll_bhhh,theta_start_pml,options);
	if twostep == 0;
	    theta_start = [theta_start_pml ; p];
		[theta_hat,FVAL,EXITFLAG,OUTPUT2] = fminunc(ll_bhhh,theta_start,options);
	end
	timetoestimate = toc(outsidetimer);

	[f,g,h] = ll_bhhh(theta_hat);
	cov = inv(h*samplesize);
    nfxp_results.title = 'Nested Fixed Point Algorithm (NFXP)';
    nfxp_results.RC(i_mc) = theta_hat(1);
    nfxp_results.c(i_mc) = theta_hat(2);
    nfxp_results.seRC(i_mc) = sqrt(cov(1,1));
    nfxp_results.sec(i_mc) = sqrt(cov(2,2));
    nfxp_results.MajorIter(i_mc) = OUTPUT1.iterations;
    nfxp_results.FuncCount(i_mc) = OUTPUT1.funcCount;
	if twostep == 0
		nfxp_results.MajorIter(i_mc) = OUTPUT1.iterations+OUTPUT2.iterations;
		nfxp_results.FuncCount(i_mc) = OUTPUT1.funcCount+OUTPUT2.funcCount;
	end
	nfxp_results.runtime(i_mc)   = timetoestimate;
	nfxp_results.converged(i_mc) = (EXITFLAG >= 1 &  EXITFLAG <= 3);
	nfxp_results.llval(i_mc)     = -FVAL*samplesize;


	% ************************************
 	% STEP 2b: ESTIMATE parameters using MPEC
 	% ************************************

	nc = 2; 				% number of cost function parameters to be estimated
	np = (twostep==0)*numel(p); 	% number of transition matrix parameters to be estimated
	[J_pattern, H_pattern] = mpec.sparsity_pattern(nc,np,mp.n, numel(p)+1);
	subplot(1,2,1), spy(J_pattern); 
	title('Jacobian of constraints')
	subplot(1,2,2), spy(H_pattern);
	title('Hessian of likelihood');

	% starting values 
	theta_start_mpec = [theta_start_pml; zeros(mp.n,1)];

	% Upper and lower bounds on parameters
	lb = zeros(2+mp.n,1); ub = zeros(2+mp.n,1);
	lb(1) = -0.000000000001; ub(1) = inf; %thetaCost
	lb(2) = -0.000000000001; ub(2) = inf; %RC

	%Put bound on EV; this should not bind, but is a cautionary step to help keep algorithm within bounds
	lb(3:end) = -5000; ub(3:end) = 0; %EV

	% No linear equality constraints
	Aeq = []; beq = [];

 	% Create mp.n dummy variables for observed mileage being equal to a particular value in the grid. 
 	% Needed for derivatives of likelihood
 	for i = 1:mp.n;
    data.xd(:,i) = (data.x == i);           
	end
	data.xd = sparse(data.xd);  % Utilize sparsity pattern. Matrix is N*T by mp.n and has only N*T nonzero elements. 

	% Define objective functions and constraints
	ll_p_mpec = @(theta) mpec.ll(data, mp, nfxp_options, theta); % objective function
	con_p_bellman = @(theta) mpec.con_bellman(data, mp, P, nfxp_options, theta); % Costraint (Bellman equation) 	

  options_mpec = optimset('DerivativeCheck','off','Display','off',...
        'GradConstr','on','GradObj','on','TolCon',1E-6,'TolFun',1E-10,'TolX',1E-15,'JacobPattern',J_pattern, 'HessPattern', H_pattern,'MaxFunEval', 100000, 'Algorithm','interior-point' ); 
  %options_mpec = optimset('GradConstr','off','GradObj','on','TolCon',1E-6,'TolFun',1E-10,'TolX',1E-15,'MaxFunEval', 100000, 'Algorithm','interior-point' ); 
  outsidetimer = tic;
  [theta_hat,FVAL,EXITFLAG,OUTPUT2] = fmincon(ll_p_mpec,theta_start_mpec,[],[],Aeq,beq,lb,ub,con_p_bellman,options_mpec);  
	timetoestimate = toc(outsidetimer);

 	mpec_results.title = 'Mathematical Programming with Equilibrium Constraints (MPEC)';
 	mpec_results.RC(i_mc) = theta_hat(1);
	mpec_results.c(i_mc) = theta_hat(2);
 	mpec_results.seRC(i_mc) = nan(1,1);
	mpec_results.sec(i_mc) = nan(1,1);
	mpec_results.MajorIter(i_mc) = OUTPUT1.iterations;
	mpec_results.FuncCount(i_mc) = OUTPUT1.funcCount;
	mpec_results.runtime(i_mc)   = timetoestimate;
	mpec_results.converged(i_mc) = (EXITFLAG == 1);
	mpec_results.llval(i_mc)     = -ll_p_mpec(theta_hat)*samplesize;

end  % End Monte Carlo

results = {nfxp_results,mpec_results};

fprintf('*************************************************************************\n'); 
fprintf('Fixed parameters and spaces\n'); 
fprintf('*************************************************************************\n'); 
fprintf('Beta           = %10.5f \n',mp.beta);
if MonteCarlo
fprintf('RC             = %10.5f \n',mp.RC);
fprintf('c              = %10.5f \n',mp.c);
end
fprintf('n              = %10.5f \n',mp.n);
fprintf('Sample size    = %10.5f \n',samplesize);
for i = 1:numel(results);
	result = results{i};
    fprintf('\n'); 
    fprintf('*************************************************************************\n'); 
    fprintf('Method, %s\n', result.title); 
    fprintf('*************************************************************************\n'); 
    fprintf('\n'); 
    fprintf('Parameter Estimates\n');
    fprintf('---------------------------------------------------------------------\n');
    fprintf('%14s %14s  %14s %14s\n', '', 'RC', 'C' ,'Likelihood');
    fprintf('---------------------------------------------------------------------\n');

    fprintf('%14s %14.3f %14.3f %14.1f \n', 'param', ...
    mean(result.RC(result.converged == 1)) , ...
    mean(result.c(result.converged == 1)) , ...
    mean(result.llval(result.converged == 1)) ); 
    fprintf('%14s %14.3f %14.3f\n\n', 's.e.', ...
    mean(result.seRC(result.converged == 1)) , ...
    mean(result.sec(result.converged == 1)) ); 
    fprintf('---------------------------------------------------------------------\n');
    fprintf('\n\n');


    fprintf('Numerical performance\n');
    fprintf('---------------------------------------------------------------------------------------\n');
    fprintf('%14s %14s  %14s %14s %14s\n', '', 'Runs Converged', 'CPU Time' ,'# of Major','# of Func.');
    fprintf('%14s %14s  %14s %14s %14s\n', '', sprintf(' (out of %g)',nMC),'(in sec.)','Iter','Eval.');
    fprintf('---------------------------------------------------------------------------------------\n');

    fprintf('%14s %14.3f %14.3f %14.1f %14.1f \n\n', '', ...
    sum(result.converged == 1) , ...
    mean(result.runtime(result.converged == 1)) , ...
    mean(result.MajorIter(result.converged == 1)), ...
    mean(result.FuncCount)); 
    fprintf('---------------------------------------------------------------------------------------\n');
end
if twostep ==0
	warning('MPEC not implemented for full MLE, only PMLE');
end
