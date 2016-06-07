function [fval,bestvalue,fminval] = wdosc20(num,numlvl,E)
 
format long g;


% User defined WDO parameters:
param.popsize = num;		% population size.
param.npar = numlvl*3;			% Dimension of the problem.
param.maxit = 700;		% Maximum number of iterations.
param.RT = 3;			% RT coefficient.
param.g = 0.2;			% gravitational constant.
param.alp = 0.4;		% constants in the update eq.
param.c = 0.4;			% coriolis effect.
maxV = 0.3;			% maximum allowed speed.
dimMin =  zeros(1,numlvl);			% Lower dimension boundary.
dimMax= 255.*ones(1,numlvl);			% Upper dimension boundary.
%---------------------------------------------------------------

% Initialize WDO population, position and velocity:
% Randomize population in the range of [-1, 1]:
pos = 2*(rand(param.popsize,param.npar)-0.5);
% Randomize velocity:
vel = maxV * 2 * (rand(param.popsize,param.npar)-0.5);  
	
%---------------------------------------------------------------

% Evaluate initial population: (Sphere Function)
for K=1:param.popsize,
	x(K,:) = (dimMax - dimMin).*((pos(K,:)+1)./2) + dimMin;
    x(K,:) = sort(round(x(K,:)),2)
    	pres(K,1) = fitnessfuncsc208(x(K,:),E);
end
%----------------------------------------------------------------

% Finding best air parcel in the initial population :
[globalpres,indx] = min(pres);
globalpos = pos(indx,:);
minpres(1) = min(pres);			% minimum pressure
%-----------------------------------------------------------------

% Rank the air parcels:
[sorted_pres rank_ind] = sort(pres);
% Sort the air parcels:
pos = pos(rank_ind,:);
keepglob(1) = globalpres;
%-----------------------------------------------------------------

% Start iterations :
iter = 1;   % iteration counter
for ij = 2:param.maxit,
    	% Update the velocity:
    	for i=1:param.popsize
		% choose random dimensions:
		a = randperm(param.npar);        			
		% choose velocity based on random dimension:
    		velot(i,:) = vel(i,a);				
        	vel(i,:) = (1-param.alp)*vel(i,:)-(param.g*pos(i,:))+ ...
				    abs(1-1/i)*((globalpos-pos(i,:)).*param.RT)+ ...
				    (param.c*velot(i,:)/i);
    	end
    
        	% Check velocity:
        	vel = min(vel, maxV);
        	vel = max(vel, -maxV);
		% Update air parcel positions:
    		pos = pos + vel;
        	pos = min(pos, 1.0);
        	pos = max(pos, -1.0); 
		% Evaluate population: (Pressure)
		for K=1:param.popsize,
			x(K,:) = sort(round((dimMax - dimMin).*((pos(K,:)+1)./2) + dimMin),2)
    			pres(K,1) = fitnessfuncsc208(x(K,:),E);
		end

    	%----------------------------------------------------
    	% Finding best particle in population
    	[minpres,indx] = min(pres);
    	minpos = pos(indx,:);           	% min location for this iteration
    	%----------------------------------------------------
    	% Rank the air parcels:
    	[sorted_pres rank_ind] = sort(pres);
    	% Sort the air parcels position, velocity and pressure:
    	pos = pos(rank_ind,:);
    	vel = vel(rank_ind,:);
    	pres = sorted_pres;  
    
    	% Updating the global best:
    	better = minpres < globalpres;
    	if better
        		globalpres = minpres;             % initialize global minimum
        		globalpos = minpos;
   	end
	% Keep a record of the progress:
    	keepglob(ij) = globalpres;
        
        minval = rank_ind(1);
%     bnvalues(ij-1,:) = x(minval,:);
fminval(ij) = globalpres;
end
minval = rank_ind(1);
bestvalue = x(minval,:)
fval = globalpres;

