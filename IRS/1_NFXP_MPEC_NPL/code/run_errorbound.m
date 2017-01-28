clear;
close all
global ev0;  ev0=0;

% Spaces
mp.n=90;				% Number of gridpoints		
mp.max=450;				% Max of mileage
mp.grid= (0:mp.n-1)';  	% Grid over mileage

% Structural parameters
mp.p=[0.0937 0.4475 0.4459 0.0127]';   	% Transition probabiliuties
mp.RC=11.7257;     						% Replacement cost
mp.c=2.45569;									% Cost parameter
mp.beta=0.9999;								% Discount factor

P0 = nfxp.statetransition(mp.p, mp.n);	% Transition matrix for mileage
cost0=0.001*mp.c*mp.grid;				% Cost function

do_nk=0;
legends={};

if do_nk
	nsaiter=5;
else
	nsaiter=10000;
end

colorOrder = get(gca, 'ColorOrder');
nfxp_options.min_cstp=nsaiter;
nfxp_options.max_cstp=nfxp_options.min_cstp;
nfxp_options.printfxp   = 0;            	% Print iteration info for fixed point algorithm if > 0
betavec=[0.95 0.99 0.999 0.9999];
for ibeta=1:numel(betavec)
	mp.beta=betavec(ibeta);
  legends{ibeta}=sprintf('beta=%g', betavec(ibeta)) ;

	[mp.ev, pk0, F, iterinfo(ibeta)]=nfxp.solve(0, P0, cost0, mp, nfxp_options); 	% Solve the model

	figure(1)
	hold on
	plot(1:numel(iterinfo(ibeta).tol), (iterinfo(ibeta).tol), 'Color', colorOrder(ibeta,:), 'LineWidth', 1.5);

	figure(2)
	hold on
	plot(1:numel(iterinfo(ibeta).tol), log(iterinfo(ibeta).tol), 'Color', colorOrder(ibeta,:),'LineWidth', 1.5);
	grid
end; 

figure(1)
xlabel('Iteration count');
ylabel('Error bound');
legend(legends, 'Location', 'SouthEast')
title('Error bound vs iteration count');
xlim([0 nsaiter-1+5*do_nk])

figure(2)
xlabel('Iteration count');
ylabel('log(Error bound)');
legend(legends, 'Location', 'SouthEast')
title('Error bound vs iteration count');
xlim([0 nsaiter-1+5*do_nk])
