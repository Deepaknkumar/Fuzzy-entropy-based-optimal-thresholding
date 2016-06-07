

function [minvalue,xmin,fminval] = BFOsc20(num,nub,E)
% ------- initialisation ----------%
Ne=10;
Nr=10;
Nc=10;
Np=num;
Ns=num/2;
D=nub*3;
C=0.01;
Ped=0.9; % elimination dispersion probability

% x=(rand(Np,D)-0.5)*60; % x lies in [-30 30]

Lb = zeros(1,D);        %%% Lower bounds
Ub = 255.*ones(1,D);    %%% Upper bounds

for i=1:Np
    xn(i,:) = Lb + (Ub-Lb).*rand(size(Lb));
end

x = sort(round(xn),2);

J=zeros(Np,1);

for k=1:Np
      J(k)=fitnessfuncsc208(x(k,:),E); % initial fitness calculation    
end

Jlast=J;

xmin = zeros(1,D);

for l=1:Ne
    for k=1:Nr
        Jchem=J;
        for j=1:Nc
            % Chemotaxis Loop %
            
            for i=1:Np
                del=(rand(1,D)-0.5)*2;
                x(i,:)=x(i,:)+(C./sqrt(del*del'))*del;
                x(i,:) = sort(round(x(i,:)),2);
                J(i)=fitnessfuncsc208(x(i,:),E); %fitness calculation
              
                for m=1:Ns
                    if J(i)<Jlast(i)
                        Jlast(i)=J(i);
                        x(i,:)=x(i,:)+C*(del./sqrt(del*del'));
                        x(i,:) = sort(round(x(i,:)),2);
                        J(i)=fitnessfuncsc208(x(i,:),E); 
                    else
                        del=(rand(1,D)-0.5)*2;
                        x(i,:)=x(i,:)+C*(del./sqrt(del*del'));
                        x(i,:) = sort(round(x(i,:)),2);
                         J(i)=fitnessfuncsc208(x(i,:),E); 
                    end   
                end
                
            end
            
            Jchem=[Jchem J];
        end  % End of Chemotaxis %
        
        
            for i=1:Np
                Jhealth(i)=sum(Jchem(i,:)); % sum of cost function of all chemotactic loops for a given k & l
            end
            
            [Jhealth1,I]=sort(Jhealth,'ascend');
            x=[x(I(1:(Np/2)),:);x(I(1:(Np/2)),:)];
            J=[J(I(1:(Np/2)),:);J(I(1:(Np/2)),:)];
%             x=[x(round(I(1:round(Np/2))),:);x(round(I(1:round(Np/2))),:)];
%             J=[J(round(I(1:round(Np/2))),:);J(round(I(1:round(Np/2))),:)];
            xmin=x(I(1),:);
    end
    Jmin(l)=min(J);
    % random elimination dispersion
   
  for i=1:Np
        r=rand;
        if r>=Ped
            x(i,:)=(rand(1,D)-0.5);
            x(i,:) = sort(round(x(i,:)),2);
            J(i)=fitnessfuncsc208(x(i,:),E);
        end
    end
       fminval(l) = Jmin(l);
end
xmin
minvalue = min(Jmin);

                