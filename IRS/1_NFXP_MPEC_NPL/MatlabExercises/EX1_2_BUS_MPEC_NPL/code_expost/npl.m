% NPL class: Solves Rust's engine repplacement model Rust(Ecta, 1987) 

classdef npl
    methods (Static) 
        function [f,g,h ]=clogit(data, v, mp, theta)
            %% Procedure for calculating the likelihood value, gradients and hessian approximation
            %  based on AM (2002) and Rust (1987).
            
            vz=v(:,1:2); % continuation from both z

            zK =[zeros(size(v,1),1),-0.001*mp.grid]+vz;
            zR =[-1,-0.001*mp.grid(1)]+vz(1,:);
            X = bsxfun(@minus,zR,zK);  
            X_data=X(data.x,:);

            [pk]=npl.pk(mp, theta, v);
            
            lp=pk(data.x);
            dk=(data.d==0);
            dr=(data.d==1);

            logl=log(lp.*dk+(1-lp).*dr);
            
            f=-mean(logl);

            NT=size(logl, 1);
        
            if nargout >=2; % compute gradient      

                res=dk-lp;                   % model "residual"
                score=zeros(NT,2);
                score(:,1)=res.*X_data(:,1); % RC
                score(:,2)=res.*X_data(:,2); % c
                g=mean(score);   

            end

             if nargout >=3; % Compute hessian
                h=score'*score/NT;
            end
        end % end of ll

        function [pk_new]=Psi(mp, pk_old, P, theta)
            %% Psi is the mapping taking old choice probabilities pk_old given current model parameters
            %  and markov transitions, and returning the optimal ccps pk_new given pk_old and parameters. 
            V      = npl.phi(mp, pk_old, P); % phi(P)
            pk_new = npl.pk(mp, theta, V);   % Psi(phi(P))
        end % end of polcy

        function [V]=phi(mp, pk, P)
            %% The function phi does most of the work for calculating phi, but it is assuming a multiplicative model
            %  so actually it only calculates
            %  V_z=(I-betaFu(P))^-1 * sum  P(a) * u(a), and
            %  V_e=(I-betaFu(P))^-1 * sum  P(a) * e(a,P), where
            %  P is actually what we call pk in the code.
            
            eulerc=0.5772156649015328606065120900824024310421;

            K_z = -mp.grid*0.001; % Maintaince cost variables (without theta multiplied!)
            R_z = -ones(size(mp.grid))+mp.grid(1)*0.001; % Replacement cost variables (without RC multiplied, so just ones!)

            K_e = eulerc - log(pk);
            R_e = eulerc - log(1-pk);

            Fu = bsxfun(@times, P, pk)  + bsxfun(@times, P(1,:), 1-pk);
            pv_z = [bsxfun(@times, R_z, 1-pk), bsxfun(@times, K_z, pk)];
            pv_e = bsxfun(@times, K_e, pk) + bsxfun(@times, R_e, 1-pk);
            
            V_z = mp.beta*P*((speye(mp.n) - mp.beta*Fu)\pv_z);
            V_e = mp.beta*P*((speye(mp.n) - mp.beta*Fu)\pv_e);
            V = [V_z V_e];

        end % end of phi

        function [pk]=pk(mp, theta, v)
            %% Calculates the choice probabilities at the grid points. 
            
            if numel(v)==1
                v=zeros(size(mp.grid,1),3);
            end
            
            v=v(:,1:2)*theta + v(:,3); % continuation form both z and e
   
            c = theta(2)*0.001*mp.grid;

            VK = - c + v;
            VR = - theta(1) - c(1) + v(1);
            pk = 1./(1+exp((VR-VK)));   
        end
           
        function [pk_solve] = solve(mp, P)         
            v   = zeros(size(mp.grid,1),3);
            pk0 = ones(mp.n,1)*0.01;
            tol = 1;
            for i_Psi = 1:100
                pk1 = npl.Psi(mp, pk0, P, [mp.RC; mp.c]);
                tol = max(abs(pk1-pk0));
                tol_vec(i_Psi) = tol;
                pk0 = pk1;
                if tol < 1e-12
                    break
                end
            end
            pk_solve = pk0;
        end
    end % end of methods
end % end of estim class