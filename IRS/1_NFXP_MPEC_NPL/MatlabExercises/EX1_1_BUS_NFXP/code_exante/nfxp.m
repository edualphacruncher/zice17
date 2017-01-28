% NFXP class: Solves Rust's engine repplacement model Rust(Ecta, 1987) 
% By Fedor Iskhakov, Bertel Schjerning, and John Rust

classdef nfxp
    methods (Static)
        function opt=setoptions(options)
            % default options
            opt.max_fxpiter= 5;             % Maximum number of times to switch between Newton-Kantorovich iterations and contraction iterations.
            opt.min_cstp   = 4;             % Minimum number of contraction steps
            opt.max_cstp   = 20;            % Maximum number of contraction steps
            opt.ctol       = 0.01;          % Tolerance before switching to N-K algorithm
            opt.rtol       = 0.02;          % Relative tolerance before switching to N-K algorithm
            opt.nstep      = 20;            % Maximum number of Newton-Kantorovich steps
            opt.ltol0      = 1.0e-10;       % Final exit tolerance in fixed point algorithm, measured in units of numerical precision
            opt.printfxp   = 1;             % print iteration info for fixed point algorithm if > 0
            opt.rtolnk     = .5;            % Tolerance for discarding N-K iteration and move to SA (tolm1/tolm > 1+opt.rtolnk)
            opt.hessian    = 'bhhh';        % Use outer product of gradient as hessian approximation

            if nargin>0
                pfields=fieldnames(options);
                for i=1:numel(pfields);
                    opt.(pfields{i})=options.(pfields{i});
                end
            end
            return
        end % end of NFXP.setoptions

        function [ev1, pk]=bellman(ev, P, c, mp)
            % NFXP.BELMANN:     Procedure to compute bellman equation
            %
            % Inputs
            %  ev0      mp.n x 1 matrix of expected values given initial guess on value function
            %
            % Outputs:
            %  ev1      mp.n x 1 matrix of expected values given initial guess of ev
            %  pk       mp.n x 1 matrix of choice probabilites (Probability of keep)
            VK=-c+mp.beta*ev;               % Value off keep
            VR=-mp.RC-c(1)+mp.beta*ev(1);   % Value of replacing    
            maxV=max(VK, VR);

            % Need to recentered by Bellman by subtracting max(VK, VR)
            ev1=P*(maxV + log(exp(VK-maxV)  +  exp(VR-maxV))); 
            
            % If requested, also compute choice probability from ev (initial input)
            if nargout>1 
                pk=1./(1+exp((VR-VK)));
            end
        end % end of NFXP.bellman                
        
        function dev=dbellman(pk,  P, mp)
            % NFXP.DBELMANN:     Procedure to compute Frechet derivative of Bellman operator
            tmp=P(:,2:mp.n).*repmat(pk(2:mp.n,1)',mp.n,1);
            dev=(mp.beta*[1-(sum(tmp,2)) tmp]);
        end % end of NFXP.dbellman
        
        function P = statetransition(p, n)
            p=[p; (1-sum(p))];
            P=0;
            for i=0:numel(p)-1;
                %keyboard
                P=P+sparse(1:n-i,1+i:n,ones(1,n-i)*p(i+1), n,n);
                P(n-i,n)=1-sum(p(1:i));
            end
            %keyboard
            P=sparse(P);
            
        end % end of NFXP.
                
        function [ev, pk, F]=solve(ev, P, cost, mp, options)
            % NFXP.SOLVE: Contraction mapping fixed point poly-algorithm.
            %
            %  syntax:	[ Pk,EV,dEV]=nfxp.solve(ev, P, cost, mp, options):
            %
            %  INPUT:
            %     ev :      m x 1 matrix or scalar zero. Initial guess of choice specific expected value function, EV.
            %     P:        State transition matrix
            %     cost:     State transition matrix
            %     mp:       Structure containing model parameters: mp.beta, mp.RC, mp.c, mp.ev
            %     options:  Optional setting for Fixed Point Algorithm
            %
            %  OUTPUT:
            %     Pk:     m x 1 matrix of conditional choice probabilities, P(d=keep|x)
            %     ev:     m x 1 matrix of expected value functions, EV(x)
            %     F:      m x m matrix of derivatives of Identity matrix I minus 
            %             contraction mapping operator, I-T' where T' refers to derivative of the expected  value function
            %
            %  FOR OPTIONS: see nfxp.setoptions
            %
            % See also:
            %   nfxp.bellman
            
            %-----------------------------------------------------------------------------------------------------------------------------
            % Default settings for Fixed point algorithm
            global BellmanIter NKIter;
            opt=nfxp.setoptions(options);
            n=mp.n;
            
            %Initialize counters
            NKIter=0;
            BellmanIter=0;
            converged=false;%initialize convergence indicator
            tolm=1;

            solutiontime=tic;

            for k=1:opt.max_fxpiter; %poly-algorithm loop (switching between SA and N-K and back)

                % SECTION A: CONTRACTION ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');
                    fprintf('Begin contraction iterations (for the %d. time)\n',k);
                    fprintf('   j           tol        tol(j)/tol(j-1) \n');
                end;
                
                %SA contraction steps
                for j=1:opt.max_cstp;
                    ev1=nfxp.bellman(ev, P, cost, mp);

                    BellmanIter=BellmanIter+1;
                    
                    tolm1=max(abs(ev-ev1));
                    rtolm=tolm1/tolm;
                    
                    if opt.printfxp>0
                        fprintf(' %3.0f   %16.8f %16.8f\n',j, tolm1,rtolm);
                    end

                    %prepare for next iteration
                    ev=ev1;
                    tolm=tolm1;

                    %stopping criteria
                    if (j>=opt.min_cstp) && (tolm1 < opt.ctol)
                        %go to NK iterations due to absolute tolerance
                        break;
                    end;
                    if (j>=opt.min_cstp) && (abs(mp.beta-rtolm) < opt.rtol)
                        %go to NK iterations due to relative tolerance
                        break
                    end
                end
                if opt.printfxp>0
                    fprintf('Elapsed time: %3.5f (seconds)\n',toc(solutiontime));
                end;
                %ev is produced after contraction steps
                %tolm will also be used below
            
                % SECTION 2: NEWTON-KANTOROVICH ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');
                    fprintf('Begin Newton-Kantorovich iterations (for the %d. time)\n',k);
                    fprintf('  nwt          tol   \n');
                    
                end
                if opt.printfxp>1
                    %plots
                    hold on
                    subplot(2,1,2), plot(ev -ev(1), '--k');
                end

                %do initial contraction iteration which is part of first N-K iteration
                [ev1, pk]=nfxp.bellman(ev, P, cost, mp); %also return choice probs=function of ev

                for nwt=1:opt.nstep; %do at most nstep of N-K steps

                    NKIter=NKIter+1;

                    %Do N-K step
                    F=speye(n) - nfxp.dbellman(pk, P,  mp);%using pk from last call to nfxp.bellman
                    ev=ev-F\(ev-ev1); %reusing ev here
                    
                    %Contraction step for the next N-K iteration
                    [ev2, pk]=nfxp.bellman(ev, P, cost, mp); %also return choice probs=function of ev

                    %Measure the tolerance
                    tolm1=max(max(abs(ev-ev2)));

                    if opt.printfxp>0
                        fprintf('   %d     %16.8e  \n',nwt,tolm1);
                    end
                    if opt.printfxp>1
                        %plot ev
                        hold on
                        subplot(2,1,2), plot(ev -ev(1), '-r');
                    end

                    %Discard the N-K step if tolm1 got worse AND we can still switch back to SA
                    if (tolm1/tolm > 1+opt.rtolnk) && (k<opt.max_fxpiter);
                        if opt.printfxp>0
                            %Discrading the N-K step
                            fprintf('Discarding N-K step\n');
                        end
                        ev=ev1; %new contraction step should start from ev1
                        break;
                    else
                        ev1=ev2; %accept N-K step and prepare for new iteration
                    end;

                    ltol=opt.ltol0;

                    if (tolm1 < ltol);
                        %N-K converged 
                        converged=true;
                        break
                    end

                end %Next N-K iteration
                if opt.printfxp>0
                    fprintf('Elapsed time: %3.5f (seconds)\n',toc(solutiontime));
                end
                if converged
                    if opt.printfxp>0
                        fprintf('Convergence achieved!\n\n');
                        % fprintf('Elapsed time: %3.5f (seconds)\n',toc(solutiontime));
                    end
                    break; %out of poly-algorithm loop
                else
                    if nwt>=opt.nstep
                        warning('No convergence! Maximum number of iterations exceeded without convergence!');
                        break; %out of poly-algorithm loop with no convergence
                    end
                end
            end
        end % end of NFXP.solve  

        function [f,g,h ]=ll(data, mp, P, options, theta)
            global ev0;

            if nargin==5;
                mp.RC=theta(1);
                mp.c=theta(2);
                if numel(theta)>2
                    mp.p=theta(3:end);
                end
            end

            mp.p=abs(mp.p); %helps BHHH which is run as unconstrained optimization
            n_c=numel(mp.c);
            N=size(data.x,1);

            % Update u, du and P evaluated in grid points
            dc=0.001*mp.grid;
            cost=mp.c*0.001*mp.grid;
            if numel(theta)>2
                P = nfxp.statetransition(mp.p, mp.n);
            end

            % Solve model
            %ev0=0;
            [ev0, pk, F]=nfxp.solve(ev0, P, cost, mp, options);

            % Evaluate likelihood function
            % log likelihood regarding replacement choice
            lp=pk(data.x);
            logl=log(lp+(1-2*lp).*(data.d));

            % add on log like for mileage process
            if numel(theta)>2
                p=[mp.p; 1-sum(mp.p)];
                n_p=numel(p)-1;
                logl=logl + log(p(1+ data.dx1));
            else
                n_p=0;
            end

            % Objective function (negative mean log likleihood)
            f=mean(-logl);

            if nargout >=2;  %% compute scores
                % step 1: compute derivative of contraction operator wrt. parameters
                dtdmp=zeros(mp.n,1+ n_c + n_p);
                dtdmp(:, 1)=P*pk-1;
                dtdmp(:, 2:1+n_c)=-(P*dc).*pk; %
                if numel(theta)>2
                    vk=-cost+mp.beta*ev0;
                    vr=-mp.RC-cost(1)+mp.beta*ev0(1);
                    vmax=max(vk,vr);
                    dtp=vmax+log(exp(vk-vmax)+exp(vr-vmax));

                    for iP=1:n_p;
                        dtdmp(1:mp.n-iP, 1+n_c+iP)= [dtp(iP:end-1)] - [dtp(n_p+1:mp.n); repmat(dtp(end), n_p-iP, 1)];
                    end
                    invp=exp(-log(p));
                    invp=[sparse(1:n_p,1:n_p,invp(1:n_p),n_p,n_p); -ones(1,n_p)*invp(n_p+1)];
                    N=size(data.x,1);
                end
                
                % step 2: compute derivative of ev wrt. parameters
                devdmp=F\dtdmp;  
                
                % step 3: compute derivative of log-likelihood wrt. parameters
                score=bsxfun(@times, (lp- 1 + data.d),[-ones(N,1) dc(data.x,:) zeros(N,n_p)] + (devdmp(ones(N,1),:)-devdmp(data.x,:)));    
                if numel(theta)>2
                    for iP=1:n_p;
                        score(:,1+n_c+iP)= score(:,1+n_c+iP) + invp(1+data.dx1,iP);
                    end
                end
                g=mean(-score,1);
            end

            if nargout >=3;  %% compute Hessian
                if strcmp(options.hessian, 'bhhh')
                    h=score'*score/(size(logl, 1)); 
                    return
                end 

                %Hessian is computed FOR STRUCTURAL PARAMETERS ONLY
                if numel(theta)>2 || n_c>1
                    error 'RTFM!'
                end

                ddiff=zeros(mp.n,2);
                % Parameters are: {RC, c}
                % Derivative of the diff between V(replace) and V(keep)
                % d(diff)/dRC
                ddiff(:,1)=-1+mp.beta*devdmp(1,1)-mp.beta*devdmp(:,1);
                % d(diff)/dc
                ddiff(:,2)=-dc(1)+mp.beta*devdmp(1,2)+dc-mp.beta*devdmp(:,2);

                % Derivative of pk
                % d(diff)/dRC
                dpk(:,1)=pk.*(pk-1).*ddiff(:,1);
                % d(diff)/dc
                dpk(:,2)=pk.*(pk-1).*ddiff(:,2);

                % Second derivatives of contraction mapping wrt parameters, RC and c
                d2t_dpi_dpj=nan(mp.n,2,2);
                % (1,1) RC RC
                d2t_dpi_dpj(:,1,1) = P*dpk(:,1);
                % (1,2) RC c
                d2t_dpi_dpj(:,1,2) = P*dpk(:,2);
                % (2,1) c RC (assuming second deriv of c =0)
                d2t_dpi_dpj(:,2,1) = -P*( dc.*dpk(:,1) + dc(1)*(1-dpk(:,1)) );
                % (2,2) c c (assuming second deriv of c =0)
                d2t_dpi_dpj(:,2,2) = -P*( dc.*dpk(:,2) + dc(1)*(1-dpk(:,2)) );


                % Cross derivative of contraction mapping wrt EV and parameters RC and c
                d2tev_dev_dpi=nan(mp.n,mp.n,2);
                % wrt RC
                tmp=P(:,2:mp.n).*repmat(dpk(2:mp.n,1)',mp.n,1);
                d2tev_dev_dpi(:,:,1)=mp.beta*[-sum(tmp,2) tmp];
                % wrt c
                tmp=P(:,2:mp.n).*repmat(dpk(2:mp.n,2)',mp.n,1);
                d2tev_dev_dpi(:,:,2)=mp.beta*[-sum(tmp,2) tmp];

                % Second derivative of contracrtion fixed point, d2ev_dpi_dpj
                % Step 1: Computed components
                q=nan(mp.n,2,2);
                for i=1:2
                    for j=1:2
                        q(:,i,j)=d2t_dpi_dpj(:,i,j)+d2tev_dev_dpi(:,:,i)*devdmp(:,j)+ ...
                                                    d2tev_dev_dpi(:,:,j)*devdmp(:,i);
                    end
                end

                d2ev_dpi_dpj=nan(mp.n,2,2);
                % Step 2: put together (separately because have to cross-reference)
                for i=1:2 %by parameter
                    d2ev_dpi_dpj(:,:,i) = F \ q(:,:,i);
                end

                H=nan(2,2);
                pkpr=pk(data.x).*(pk(data.x)-1);
                H(1,1)=sum(pkpr.*ddiff(data.x,1).*ddiff(data.x,1) + mp.beta*(lp-1+data.d).*(d2ev_dpi_dpj(1,1,1)-d2ev_dpi_dpj(data.x,1,1)));
                H(2,2)=sum(pkpr.*ddiff(data.x,2).*ddiff(data.x,2) + mp.beta*(lp-1+data.d).*(d2ev_dpi_dpj(1,2,2)-d2ev_dpi_dpj(data.x,2,2)));
                % H(1,2)=sum(pkpr.*ddiff(data.x,1).*ddiff(data.x,2) + mp.beta*(lp-1+data.d).*(d2ev_dpi_dpj(1,1,2)-d2ev_dpi_dpj(data.x,1,2)));
                H(1,2)=sum(pkpr.*ddiff(data.x,2).*ddiff(data.x,1) + mp.beta*(lp-1+data.d).*(d2ev_dpi_dpj(1,2,1)-d2ev_dpi_dpj(data.x,2,1)));
                H(2,1)=H(1,2);

                h=-H/numel(logl);
            end
        end % end of NFXP.ll

        function [data] = simdata(N,T,mp, P, pk)
            if nargin ==3  % solve model if solution is not provided
                fxp_opt.printfxp=0;
                P = nfxp.statetransition(mp.p, mp.n);
                cost=0.001*mp.c*mp.grid;
                [ev, pk]=nfxp.solve(0, P, cost, mp, fxp_opt);
            end
            
            id=repmat((1:N),T,1);
            t=repmat((1:T)',1,N);
            u_dx=rand(T,N);
            u_d=rand(T,N);
            
            csum_p=cumsum(mp.p);
            dx1=0;
            for i=1:numel(csum_p);
                dx1=dx1+ (u_dx>(csum_p(i)));
            end;
            
            x=nan(T, N);
            x1=nan(T, N);
            x(1,:)=ones(1, N); % Intial conditions
            
            for it=1:T;
                d(it,:)=(u_d(it,:)<(1-pk(x(it,:)')'));  % Replace = 1, Keep=0
                x1(it,:)=min(x(it,:).*(1-d(it,:)) + d(it,:) + dx1(it,:), mp.n); % odometer reading next period
                if it<T;
                    x(it+1,:)=x1(it,:);
                end
            end
            
            data.id=id;
            data.t=t;
            data.d=d;
            data.x=x;
            data.x1=x1;
            data.dx1=dx1;
            
            pfields=fieldnames(data);
            for i=1:numel(pfields);
                data.(pfields{i})=reshape(data.(pfields{i}), T*N,1);
            end
        end % end of NFXP.simdata
        
        function dta = readbusdata(mp,bustypes)
            load('busdata1234.mat');

            % Select busses
            if nargin>1
                data=data(data(:,2) <=bustypes,:);
            end

            id=data(:,1);       % Bus id
            bustype=data(:,2);  % bustype: 1,2,3,4
            d1=data(:,5);       % Lagged replacement dummy, d_(t-1)
            d=[d1(2:end);0];    % Replacement dummy, d_t
            x=data(:,7);        % Odometer, x_t

            %keyboard
            % Discretize odometer data into 1, 2, ..., n
            x=ceil(x.*mp.n/(mp.max*1000));                                 

            % Mothly milage
            dx1=x(1:end,1)-[0;x(1:end-1,1)];           % Take first difference on odometer                           
            dx1=dx1.*(1-d1)+x.*d1;                     % Make true x_t-x_(t-1) (replace first difference by x_t if replacement dummy lagged is 1)

            % remove observations with missing lagged mileage
            data=[id, bustype, d,x,dx1];                                           % [ id , replace dummy lagged, odometer, "change in odometer" , bustype]
            remove_first_row_index=data(:,1)-[0;data(1:end-1,1)];                  % Take first difference of ID...
            data=data((remove_first_row_index==0),:);                              % ... and only keep lines where ID hasn't changed 

            % Save data structure
            dta.d=data(:,3);           
            dta.x=data(:,4);
            dta.dx1=data(:,5);
        end % end of NFXP.readbusdata    


      
   
    end % end of methods
    
    
end % end of nfxp class
