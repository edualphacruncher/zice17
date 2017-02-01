% setup.m: sets up globals for solution of mixed strategy equilibrium in simultaneous
%          move leapfrogging game
%          JOhn Rust, Georgetown University, Jan  2017

global bet c0 nstates cgrid vmat k k1 k2 dtp;

dtp=0;      % 1 for deterministic technological progress  (see stp.m)
c0=5;       % top apex cost state for game
nstates=25;  % number of discrete cost states 
bet=.9523;    % discount factor
k=5.2;
k1=5.2;
k2=0;

cgrid=(0:c0/(nstates-1):c0)';
%cgrid=(0:.1/(nstates-1):.1)';

vmat=[];
