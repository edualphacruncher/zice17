classdef mpec
    % MPEC class for structural estimation of discrete choice models.
    % By Bertel Schjerning & Patrick Kofod Mogensen
    methods (Static)
        function [f,g,h ]=ll(data, mp, options, theta)
            global ev0;
            mp.RC=theta(1);
            mp.c=theta(2);
            if numel(theta)>(2+mp.n)
                n_p=numel(mp.p);
                mp.p=theta(3:3+n_p)';
                ev=theta(4+n_p:end);
            else
                n_p=0;
                ev=theta(3:end);
            end

            % Update u, du and P evaluated in grid points
            dc=0.001*mp.grid;
            c=mp.c*0.001*mp.grid;

            VK=-c+mp.beta*ev;
            VR=-mp.RC-c(1)+mp.beta*ev(1);
            pk=1./(1+exp((VR-VK)));
            lp=pk(data.x);
            dk=(data.d==0);
            dr=(data.d==1);

            logl=log(dk.*lp+(1-lp).*dr);

            % add on log like for mileage process
            if n_p>0;
                logl=logl + log(mp.p(1+ data.dx1)');
            end

            f=-mean(logl);

            if nargout >=2; % compute gradient        
                res=lp-dk; % model "residual"
                g=zeros(2+mp.n,1);
                g(1)=-mean(res);                            %RC
                g(2)=mean(res.*(dc(data.x)-dc(1)));         %c
                g(3)=-mp.beta*mean(res.*(data.xd(:,1)-1));    %ev(1)

                NT=numel(res);
                for i=2:mp.n;        
                    g(2 + n_p+ i)=-mp.beta*sum(res.*(data.x==i))/NT;         %ev(2:N)
                    g(2 + n_p+ i)=-mp.beta*sum(res(data.x==i))/NT;         %ev(2:N)
                    g(2 + n_p+ i)=-mp.beta*sum(res(data.xd(:,i)))/NT;      %ev(2:N) sparse index precalculated
                end
                g=-g;
            end

            if nargout >=3; % Compute hessian
                %Hessian is computed for only RC, c, and ev
                if n_p>0
                    error 'RTFM!'
                end

                % Parameters are: {RC, c, ev}
                % Derivative of the diff between V(replace) and V(keep)
                ddiff=zeros(mp.n,2+mp.n);
                ddiff(:,1)=-1;                          % d(diff)/dRC
                ddiff(:,2)=dc-dc(1);                    % d(diff)/dc
                ddiff(2:end,3)=mp.beta*ones(mp.n-1,1);   % d(diff)/ev(1)
                ddiff(2:end,4:end)=-mp.beta*eye(mp.n-1);

                %Compute Hessian of likelihood [2x2]
                H=sparse(2+mp.n,2+mp.n);
                for i=1:2+mp.n
                     for j=i:2+mp.n;
                        if ((i<=3)||(i==j))
                             % TODO: skip elements (sparsity), utilize symmetry and precompute indices
                             H(i,j)=sum(pk(data.x).*(pk(data.x)-1).*ddiff(data.x,i).*ddiff(data.x,j));
                             H(j,i)=H(i,j);
                        end
                    end
                end       
                h=-H/NT; % since we are minimizing the negative of the mean log-likleihood
            end
         end % end of ll
         function [c,ceq,DC,DCeq] = con_bellman(data, mp, P, par, theta)
            mp.RC=theta(1);
            mp.c=theta(2);
            if numel(theta)>(2+mp.n)
                error('RTFM! con_bellman not implemented for full mle')
            else
                n_p=0;
                ev=theta(3:end);
            end

            % Update u, du and P evaluated in grid points
            dc=       0.001*mp.grid;
            cost=mp.c*0.001*mp.grid;

            [ev1, pk]=nfxp.bellman(ev, P, cost, mp);

            % Define and evaluate nonlinear inequality constraints
            c = [];

            % Define and evaluate nonlinear equality constraints
            ceq = ev1-ev;

            % Define and evaluate the constraint Jacobian (DC, DCeq).   
            if nargout >= 3
                DC= [];
                DCeq=[];

                DCeq=sparse(mp.n,2+mp.n); 
                tmp=(P*pk);
                DCeq(:,1)=- P*(1*(1-pk));
                DCeq(:,2)=-P*((dc-dc(1)).*pk);

                dev=nfxp.dbellman(pk,  P, mp);
                DCeq(:,3:end)=(dev-speye(mp.n));
                DCeq=DCeq';
            end
            
   
        end
             function [JacobSpaPattern, HessSpaPattern] = sparsity_pattern(nc,np,N, M)
             BandDiag=0;
             for i=0:M-1
                 a=diag(ones(N,1),i);
                 BandDiag=BandDiag+a(1+i:end,1+i:end);
             end
             BandDiag(:,1)=ones(N,1);     
             JacobSpaPattern=[ones(N, nc+np) BandDiag];	

             HessSpaPattern = ones(nc+np+1,nc+np+N);
             HessSpaPattern = [HessSpaPattern; ones(N-1,nc+np+1) eye(N-1)];

          end
    end
end

