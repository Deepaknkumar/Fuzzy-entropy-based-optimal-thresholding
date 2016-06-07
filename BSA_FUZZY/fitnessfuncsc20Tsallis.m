function z = fitnessfuncsc20Tsallis(u,E)
%             p = E;
            u
            threshlevel = length(u);
            q = 10;
            %u = round(u);
            omega = 0;
            %%%---------------- For First threshold level--------------%%%
            omega = sum(E(1:u(1)));         %%% Calculation of omega for first level
%             S(1) = sum( (E(1:u(1))./omega).*(log(E(1:u(1))./omega)) );
            sterm = (E(1:u(1))/omega);
%             sterm = sum(E(1:u(1))/omega);
            S(1) = (1-sum((sterm).^q))/(q-1);  
                for ij = 2:threshlevel          %%% Calculation of Entropy of each level
                   omega = sum(E(u(ij-1)+1:u(ij)));  
                   sterm = (E(u(ij-1)+1:u(ij))/omega);
%                    sterm = sum(E(u(ij-1)+1:u(ij))/omega);
                   S(ij) = (1-sum((sterm).^q))/(q-1);  
                end
            omega = sum(E(u(ij)+1:256));
            sterm = (E(u(ij)+1:256)/omega);
%             sterm = sum(E(u(ij)+1:256)/omega);
            S(ij+1) = (1-sum((sterm).^q))./(q-1); 
%             S(ij+1) = sum( (E(u(ij)+1:256)./omega).*(log(E(u(ij)+1:256)./omega)) );
            z = -1*(1*(sum(S) + (1-q)*prod(S)));       %%% Sum of entropies of all levels
            
            threshvalues = [1 u 256];
            H = [];
            cdf = [];
            prob = E;
            cdf(1) = prob(1);
            
            for i=1:255
                cdf(i+1) = cdf(i) + prob(i+1);
            end
            
            for i=1:256
                if cdf(i)==0
                    cdf(i) = 1;
                end
            end
            
            