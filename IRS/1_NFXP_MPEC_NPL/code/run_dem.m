clear; close all;
global ev0;  ev0=0;

%% NFXP options
nfxp_options=nfxp.setoptions;  			% Set default options
nfxp_options.printfxp   = 1            	% Print iteration info for fixed point algorithm if > 0

% nfxp_options.min_cstp=100;
% nfxp_options.max_cstp=nfxp_options.min_cstp;

% Spaces
mp.n=90;				% Number of gridpoints		
mp.max=450;				% Max of mileage
mp.grid= (0:mp.n-1)';  	% Grid over mileage

% Structural parameters, model 11
for dynamic=0:1;
if dynamic
	mp.RC=10.0750;
	mp.c=2.2930;
	mp.p=[0.3919 0.5953 ]'	
	mp.beta=0.9999;					
	pl='-b';				
else
	mp.RC=7.6358;
	mp.c=71.5133;
	mp.p=[0.3919 0.5953 ]'	
	mp.beta=0;
	pl='-r';				
end
mp0=mp;

P = nfxp.statetransition(mp.p, mp.n);	% Transition matrix for mileage
cost=0.001*mp.c*mp.grid;				% Cost function

[ev0, pk]=nfxp.solve(ev0, P, cost, mp, nfxp_options); 	% Solve the model
[pp, pp_K, pp_R] = nfxp.eqb(mp,P, pk); 

	x=mp.max*1000*(1:mp.n)'; 
	fprintf('Fraction of bus engines replaced each month : %1.5f \n', sum(pp_R'));
	fprintf('Mean lifetime of bus engine (months)        : %1.5f \n', 1/sum(pp_R'));
	fprintf('Mean mileage at overhaul                    : %1.5f \n',  sum(x.*pp_R')/(mp.n*sum(pp_R')));
	fprintf('Mean mileage since last Replacement         : %1.5f \n',  sum(x.*pp_K')/(mp.n*sum(pp_K')));

figure(1)
plot(mp.grid, [pp_K' 100*pp_R']);
legend('Pr(x, i=Keep)', 'Pr(x, i=Replace)')
title('Equilibrium Distribution: Bus mileage');
xlabel('Mileage');
ylabel('CDF'); 



RCgrid=1:0.5:30;
nfxp_options.printfxp=0;
for i=1:numel(RCgrid);
	mp.RC=RCgrid(i);
	[ev0, pk]=nfxp.solve(ev0, P, cost, mp, nfxp_options); 	% Solve the model
	[pp, pp_K, pp_R] = nfxp.eqb(mp,P, pk);   % compute equilibrium distribution
	Demand(i)=12*sum(pp_R');	% Demand: Expected number of bus engines replaced (in equilibrium)
end

figure(2)
hold on;
plot(RCgrid*7513/mp0.RC, Demand, pl);
ylim([0 1])
xlim([0 12000])
title('Expected Replacement Demand Function');
xlabel('Replacement cost, RC');
ylabel('Expected Annual Engine Replacement'); 
legend('beta==0', 'beta=0.9999')
end
