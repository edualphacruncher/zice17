global ev0;  ev0=0;
%% NFXP options
nfxp_options=nfxp.setoptions;  			% Set default options
nfxp_options.printfxp   = 0;            % Print iteration info for fixed point algorithm if > 0
% nfxp_options.hessian  = 'analytical';		

% Print nfxp options: 
nfxp_options=nfxp_options
nfxp_options.nstep=1

% Swiches:
MonteCarlo=1;   % 1: Do monte Carlo, 0: Use rust bus data 
twostep=0;      % 1: Estimate using two-step pmle, 0: Use full mle
bustypes=4;		  % Select bus types smaller than this number, i.e. 1,2,..bustypes

% Spaces
mp.n=175;				% Number of gridpoints		
mp.max=450;				% Max of mileage
mp.grid= (0:mp.n-1)';  	% Grid over mileage

% Structural parameters
mp.p=[0.0937 0.4475 0.4459 0.0127]';   	% Transition probabiliuties
mp.RC=11.7257;     											% Replacement cost
mp.c=2.45569;														% Cost parameter
mp.beta=0.9999;													% Discount factor

if MonteCarlo     % Solve model for DGP
	nMC=1; 	  % number of MC samples
	N=50;		  % Number of busses to simulate 
	T=119;		  % Number of time periods to simulate 

	P0 = nfxp.statetransition(mp.p, mp.n);	% Transition matrix for mileage
	cost0=0.001*mp.c*mp.grid;				% Cost function
	[mp.ev, pk0]=nfxp.solve(0, P0, cost0, mp, nfxp_options); 	% Solve the model

else    % Read in and discretize Rust's data
	nMC=1; 
 	timetosimulate=tic;
	data = nfxp.readbusdata(mp, bustypes);
	timetosimulate=toc(timetosimulate);
end

% ************************************************************************
if MonteCarlo;
	fprintf('Begin Monte Carlo, with n=%d replications\n', nMC);
else
	fprintf('Structural Estimation using data from Rust(1987)\n', nMC);
end
% ************************************************************************
rand('seed',300);

%main loop
for i_mc=1:nMC;
	% ************************************
	% STEP 0: SIMULATE DATA 
	% ************************************

	if MonteCarlo;
		timetosimulate=tic;
		data = nfxp.simdata(N, T, mp, P0, pk0);
		timetosimulate=toc(timetosimulate);
		fprintf('i_mc=%d, Time to simulate data : %1.5g seconds\n', i_mc, timetosimulate);
	end
	samplesize=numel(data.d);
	% ************************************
	% STEP 1: ESTIMATE p 
	% ************************************

	% Estimate model using Partial MLE
	tab=tabulate(data.dx1);
	tab=tab(tab(:,3)>0,:);
	p=tab(1:end-1,3)/100;
	P = nfxp.statetransition(p, mp.n);

  % ************************************
	% STEP 2: ESTIMATE parameters 
	% ************************************

	% Starting values
	mp.beta=mp.beta;   
	mp.p=p;

    theta_start=[0;0];
    ll_bhhh=@(theta)   nfxp.ll(data, mp, P, nfxp_options, theta);
    
	options= optimset('Display','off', 'GradObj','on', 'TolFun',1E-8,'TolX',0,'Hessian','on');

    outsidetimer=tic;
	[theta_hat,FVAL,EXITFLAG,OUTPUT1]=fminunc(ll_bhhh,theta_start,options);
	if twostep==0;
	    theta_start=[theta_hat ; p];
		[theta_hat,FVAL,EXITFLAG,OUTPUT2]=fminunc(ll_bhhh,theta_start,options);
	end
	timetoestimate=toc(outsidetimer);

	[f,g,h]=ll_bhhh(theta_hat);
	cov=inv(h*samplesize);
 	bhhh.RC(i_mc)=theta_hat(1);
	bhhh.c(i_mc)=theta_hat(2);
 	bhhh.seRC(i_mc)=sqrt(cov(1,1));
	bhhh.sec(i_mc)=sqrt(cov(2,2));
	bhhh.MajorIter(i_mc)=OUTPUT1.iterations;
	bhhh.FuncCount(i_mc)=OUTPUT1.funcCount;
	if twostep==0
		bhhh.MajorIter(i_mc)=OUTPUT1.iterations+OUTPUT2.iterations;
		bhhh.FuncCount(i_mc)=OUTPUT1.funcCount+OUTPUT2.funcCount;
	end
	bhhh.runtime(i_mc)=timetoestimate;
	bhhh.converged(i_mc)	= 	(EXITFLAG>=1 &  EXITFLAG<=3);
	bhhh.llval(i_mc)	= 	-FVAL*samplesize;



end  % End Monte Carlo


fprintf('\n');
fprintf('Beta           = %10.5f \n',mp.beta);
if MonteCarlo
	fprintf('RC             = %10.5f \n',mp.RC);
	fprintf('c              = %10.5f \n',mp.c);
end
fprintf('n              = %10.5f \n',mp.n);
fprintf('Sample size    = %10.5f \n',samplesize);
fprintf('\n'); 
fprintf('\n'); 
fprintf('Parameter Estimates\n');
fprintf('---------------------------------------------------------------------\n');
fprintf('%14s %14s  %14s %14s\n', '', 'RC', 'C' ,'Likelihood');
fprintf('---------------------------------------------------------------------\n');

result=bhhh;
fprintf('%14s %14.3f %14.3f %14.1f \n', 'param', ...
mean(result.RC(result.converged==1)) , ...
mean(result.c(result.converged==1)) , ...
mean(result.llval(result.converged==1)) ); 
fprintf('%14s %14.3f %14.3f\n\n', 's.e.', ...
mean(result.seRC(result.converged==1)) , ...
mean(result.sec(result.converged==1)) ); 
fprintf('---------------------------------------------------------------------\n');
fprintf('\n\n');


fprintf('Numerical performance\n');
fprintf('-----------------------------------------------------------------------------------------------------------------------\n');

fprintf('%14s %14s  %14s %14s %14s\n', '', 'Runs Converged', 'CPU Time' ,'# of Major','# of Func.');
fprintf('%14s %14s  %14s %14s %14s\n', '', sprintf(' (out of %g)',nMC),'(in sec.)','Iter','Eval.');
fprintf('-----------------------------------------------------------------------------------------------------------------------\n');

result=bhhh;
fprintf('%14s %14.3f %14.3f %14.1f %14.1f \n\n', '', ...
sum(result.converged==1) , ...
mean(result.runtime(result.converged==1)) , ...
mean(result.MajorIter(result.converged==1)), ...
mean(result.FuncCount)); 
fprintf('-----------------------------------------------------------------------------------------------------------------------\n');

