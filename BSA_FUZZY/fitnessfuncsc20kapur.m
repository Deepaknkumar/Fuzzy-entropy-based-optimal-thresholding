function z = fitnessfuncsc20kapur(u,E)
%             threshlevel = length(u);
%             %u = round(u);
%             for i=1:length(E)
%                 if E(i)==0
%                     E(i)=1;
%                 end
%             end
%             omega = 1;
%             %%%---------------- For First threshold level--------------%%%
%             omega = sum(E(1:u(1)));         %%% Calculation of omega for first level
% %             if omega==0
% %                 omega=1
% %             end
%             H(1) = sum( (E(1:u(1))./omega).*(log(abs(E(1:u(1))./omega))) );
%                 for ij = 2:threshlevel          %%% Calculation of Entropy of each level
%                    omega = sum(E(u(ij-1)+1:u(ij)));     
% %                    if omega==0
% %                        omega=1
% %                    end
%                    H(ij) = sum( (E(u(ij-1)+1:u(ij))./omega).*(log(E(u(ij-1)+1:u(ij))./omega)) ); 
%                 end
%             omega = sum(E(u(ij)+1:256));
% %             if omega==0
% %                 omega=1
% %             end
%             H(ij+1) = sum( (E(u(ij)+1:256)./omega).*(log(E(u(ij)+1:256)./omega)) );
%             z = -sum(H);         %%% Sum of entropies of all levels
            
u
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
            
            for i=1:length(threshvalues)-1
                if i==1
                    current_pdf = prob(threshvalues(i):threshvalues(i+1));
                else
                    current_pdf = prob(threshvalues(i)+1:threshvalues(i+1));
                end
                if threshvalues(i+1) == 0
                    threshvalues(i+1) = 1;
                end
                if threshvalues(i) == 0
                    threshvalues(i) = 1;
                end
                
                Psum = cdf(threshvalues(i+1)) - cdf(threshvalues(i));
                
                if Psum == 0
                    Psum = 1;
                end
                logterm = (current_pdf/Psum);
                for k=1:length(logterm)
                    if logterm(k)==0
                        logterm(k)=.00000000000000001;
                    end
                end
                H(i) = sum((logterm).*log(logterm));
            end
            
            z=sum(H);