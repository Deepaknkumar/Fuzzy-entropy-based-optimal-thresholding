function z = fitnessfuncsc208(u,prob)
            for i=1:256
                if prob(i)==0
                    prob(i) = 1;
                end
            end
            
            Nparts = (length(u)/3)+1;
            
            probmat = [];
            for i=1:Nparts
                probmat(i,:) = prob;
            end
%             mu1 = zeros(1,256);
%             mu12 = zeros(1,256);
%             mu23 = zeros(1,256);
%             mu34 = zeros(1,256);
%             mu45 = zeros(1,256);
%             mu56 = zeros(1,256);
%             mu67 = zeros(1,256);
%             mu78 = zeros(1,256);
%             mu8 = zeros(1,256);
            
%             mu = [];
            
            
%             for i=1:Nparts
%                 mu(i,:) = zeros(1,256);
%             end
%             
            mu = zeros(Nparts,256);
            
            mu(1,:) = mucreator([0 0 0 u(1:3)]);
            for i=2:Nparts-1
                mu(i,:) = mucreator(u(((i-2)*3)+1:3*i));
            end
            mu(end,:) = mucreator([u(((Nparts-2)*3)+1:end) 256 256 256]);
            
%             mu1 = mucreator([0 0 0 u(1:3)]);
%             mu12 = mucreator(u(1:6));
%             mu23 = mucreator(u(4:9));
%             mu34 = mucreator(u(7:12));
%             mu45 = mucreator(u(10:15));
%             mu56 = mucreator([u(13:15) 256 256 256]);
            

            probsums = [];
%             for i=1:Nparts
%                 probsums(i) = sum(prob.*mu(i,:));
%             end
            probsums = sum((probmat.*mu),2);
%             prob1 = sum(prob.*mu1);
%             prob12 = sum(prob.*mu12);
%             prob23 = sum(prob.*mu23);
%             prob34 = sum(prob.*mu34);
%             prob45 = sum(prob.*mu45);
%             prob5 = sum(prob.*mu56);
            
%             H1 = 0;
%             H12 = 0;
%             H23 = 0;
%             H34 = 0;
%             H45 = 0;
%             H5 = 0;
            
            H = zeros(Nparts,1);
            
            ltmu = [];          %% log term mu
            for i=1:Nparts
                ltmu(i,:) = (prob.*mu(i,:))/probsums(i);
            end
            
%             ltmu1 = (prob.*mu1)/prob1;
%             ltmu12 = (prob.*mu12)/prob12;
%             ltmu23 = (prob.*mu23)/prob23;
%             ltmu34 = (prob.*mu34)/prob34;
%             ltmu45 = (prob.*mu45)/prob45;
%             ltmu5 = (prob.*mu56)/prob5;
            
            %%% Making zeros .000000000000001
            zeromatbin = ltmu==0;
            conversionmat = (.000000000000001).*zeromatbin;
            ltmu = ltmu + conversionmat;
%             for i=1:256
%                 if ltmu1(i)==0
%                     ltmu1(i) = .00000000000000001;
%                 end
%                 if ltmu12(i)==0
%                     ltmu12(i)= .0000000000000001;
%                 end
%                 if ltmu23(i)==0
%                     ltmu23(i) = .00000000000000001;
%                 end
%                 if ltmu34(i)==0
%                     ltmu34(i) = .00000000000000001;
%                 end
%                 if ltmu45(i)==0
%                     ltmu45(i)= .0000000000000001;
%                 end
%                 if ltmu5(i)==0
%                     ltmu5(i) = .00000000000000001;
%                 end
%             end
            
%             for i=1:256
%                 H1 = H1 + (prob(i)*mu1(i)*log(ltmu1(i))/prob1);
%                 H12 = H12 + (prob(i)*mu12(i)*log(ltmu12(i))/prob12);
%                 H23 = H23 + (prob(i)*mu23(i)*log(ltmu23(i))/prob23);
%                 H34 = H34 + (prob(i)*mu34(i)*log(ltmu34(i))/prob34);
%                 H45 = H45 + (prob(i)*mu45(i)*log(ltmu45(i))/prob45);
%                 H5 = H5 + (prob(i)*mu56(i)*log(ltmu5(i))/prob5);
%             end
           
            for i=1:256
                H = H + ((probmat(:,i).*mu(:,i).*log(ltmu(:,i)))./probsums);
            end
%             for i=1:256
%                 Hd = Hd + (prob(i)*mud(i)/probD);
%                 Hm = Hm + (prob(i)*mum(i)/probM);
%                 Hb = Hb + (prob(i)*mub(i)/probB);
%             end
            zsumval = 0;
            for i=1:Nparts
                if i==Nparts
                    zsumval = zsumval + abs(H(Nparts)-H(1));
                else
                zsumval = zsumval + abs(H(i)-H(i+1));
                end
            end
%             z = abs(H1- H12) + abs(H12-H23) + abs(H23-H34) + abs(H34-H45) + abs(H45-H5) + abs(H5-H1); %+ -(abs(u(1)-0) + abs(u(2)-u(1)) + abs(u(3)-u(2)) + abs(u(4)-u(3)) + abs(u(5)-u(4)) + abs(u(6)-u(5)) + abs(255-u(6)));
            z = zsumval;