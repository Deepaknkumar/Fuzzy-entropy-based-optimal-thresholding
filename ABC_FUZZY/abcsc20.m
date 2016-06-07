function [minvalue,xo,fminval] = abcsc20(num,nub,E)
SN=num;
limit=10;
n=nub*3;
cycle=700;
Xminj = zeros(1,n);   % lower bound for dimen j
Xmaxj = 255.*ones(1,n);   % upper bound for dimen j

for x=1:SN
    for y=1:n
    X(x,y)= Xminj(1,y) + (rand(1)*(Xmaxj(1,y)-Xminj(1,y)));
    end
    X(x,:) = sort(round(X(x,:)),2);
    trial(x)=0;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% make in order
    for g=1:n-1
        min = X(x,g);
        for l=g:n
            if X(x,l) <= min
                    min = X(x,l);
                 si=l;
            end
         end
    temp=X(x,g);
    X(x,g)=X(x,si);
    X(x,si)=temp;
    
    end
end



%   Evaluate Population
for s=1:SN    
 z = fitnessfuncsc208(X(s,:),E);         %evaluate Fitness
 fX(s)=z;                           %store in a vector
end
        maxim = fX(1);
    for l=1:SN                      %find Maximum Fitness value of food
        if fX(l) >= maxim
             maxim = fX(l);
             si=l;
        end
    end
    
fo=fX(si) ;          %max fitness
xo=X(si,:) ;         %best food source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DEFINE FO AND XO

V=X;
for a=1:cycle
    
    for i=1:SN
        j=ceil(rand *(n-1));            %for Random Bee
        k=ceil(rand *(SN-1));          %for random food source
          
                V(i,j)=X(i,j)+ rand(1)*(X(k,j)-X(i,j)); %%%% for employed bee - convert some bees in employed ones
                V(i,j) = sort(round(V(i,j)),2);
                %%arrange in ascending order
                for g=1:n-1
                    min = V(i,g);
                for l=g:n
                    if V(i,l) <= min
                         min = V(i,l);
                         si=l;
                    end
                end
                
                temp=V(i,g);
                V(i,g)=V(i,si);
                V(i,si)=temp;
                end
                %%arrangement done
                
                fv(i)=fitnessfuncsc208(V(i,:),E);        %calculate Fitness of newly found food source
                
                if fv(i)< fo                    %compare fitness with current maximum
                     X(i,:)=V(i,:);             %if greater replace the original vector with this new one
                     xo=X(i,:);                  %update best(maximum) solution
                     fo=fv(i);                    %update max fitness
                 else
                     trial(i)=trial(i)+1;         %update trial
                 end
     end       
              
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    for i=1:SN              %for complete matrix of solution
        s=0;                
        for j=1:SN          %sum of fitness of all the rows
            s=s+fX(j);
        end
        P(i)=fX(i)/s;       %probablity distribution
    end
       
    t=0; 
    for i=1:SN
    for R=t:SN
        
        if rand < P(i)
                j=ceil(rand *(n-1));
                k=ceil(rand *(SN-1));
                V(i,j)= X(i,j) + rand(1)*( X(k,j)-X(i,j) ) ;     %% for onlooker
                V(i,j)= sort(round(V(i,j)),2);
                for g=1:n-1
                    min = V(i,g);
                for l=g:n
                    if V(i,l) <= min
                         min = V(i,l);
                         si=l;
                    end
                end
                
                temp=V(i,g);
                V(i,g)=V(i,si);
                V(i,si)=temp;

                end
                
                fv(i)=fitnessfuncsc208(V(i,:),E);
            
                if fv(i)< fo
                 X(i,:)=V(i,:);
                 xo=X(i,:);
                 fo=fv(i);
             else
                 trial(i)=trial(i)+1;
             end  
        end
        
            
    end
        if (trial(i))>limit
        for j=1:n
        X(i,j)= Xminj(1,j) + rand(1)*(Xmaxj(1,j)-Xminj(1,j));
        end
        X(i,:) = sort(round(X(i,:)),2);
    end
    end
   fminval(a) = fo; 
end

for g=1:n-1
        min = xo(g);
    for l=g:n
        if xo(l) <= min
             min = xo(l);
             si=l;
        end
    end
    temp=xo(g);
    xo(g)=xo(si);
    xo(si)=temp;
 end

xo
minvalue = fo;    
