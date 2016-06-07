function [fmin,bestnest,fminval] = cuckoosc20(n,num,E)

if nargin<1
    n = 25;
end

format long;

number_of_solution = num*3;
% Lb = 0.*ones(1,number_of_solution);
% [~,lowerindx] = min(E);
Lb = 0.*ones(1,number_of_solution);
Ub = 255.*ones(1,number_of_solution);

for i=1:n
    nest(i,:) = Lb + (Ub - Lb).*rand(size(Lb));
    nest(i,:) = sort(round(nest(i,:)),2);
end

pa = .25;           %discovery rate

fitness = 10^10.*ones(n,1);
[fmin,bestnest,nest,fitness] = get_best_nest(nest,nest,fitness,E);

N_iterTotal = 1000; %total number of iterations

N_iter = 0;

for iter=1:N_iterTotal
    
    new_nest = get_cuckoos(nest,bestnest,Lb,Ub);
    [fnew,best,nest,fitness] = get_best_nest(nest,new_nest,fitness,E);
    
    N_iter = N_iter + n;
    
    new_nest = empty_nests(nest,Lb,Ub,pa);      %discovery and Randomization
    [fnew,best,nest,fitness] = get_best_nest(nest,new_nest,fitness,E);
    
    N_iter = N_iter + n;
    
    if fnew < fmin
        fmin = fnew;
        bestnest = best;
    end
%     clc
%     fprintf('\n Cuckoo Search %d %% Completed......\n',uint8((iter*100)/N_iterTotal))
    fminval(iter) = fmin;
end

bestnest


%subfunctions

%cuckoos using levy flight
function nest = get_cuckoos(nest,best,Lb,Ub)
n = size(nest,1);
beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n
    s=nest(j,:);
    u =randn(size(s))*sigma;
    v =randn(size(s));
    step = u./abs(v).^(1/beta);
    
    stepsize = .01*step.*(s - best);
    
    s = s + stepsize.*randn(size(s));

    nest(j,:) = sort(round(SimpleBounds(s,Lb,Ub)),2);
end

%current best nest

function [fmin,best,nest,fitness] = get_best_nest(nest,newnest,fitness,E)

for j = 1:size(nest,1)
    fnew = fobj(newnest(j,:),E);
    if fnew < fitness(j)
        fitness(j) = fnew;
        nest(j,:) = newnest(j,:);
    end
end
%current bests
[fmin,K] = min(fitness);
best = nest(K,:);

%replacing some nests by constructing new solutions/nests
function new_nest = empty_nests(nest,Lb,Ub,pa)
n = size(nest,1);

k = rand(size(nest)) > pa;
stepsize = rand*(nest(randperm(n),:) - nest(randperm(n),:));
 new_nest = nest + stepsize.* k ;
 
 for j=1:size(new_nest)
     s = new_nest(j,:);
     new_nest(j,:) = sort(round(SimpleBounds(s,Lb,Ub)),2);
 end
 
 %function for simple bounds
    function s  = SimpleBounds(s,Lb,Ub)
        ns_temp = s;
        I = s<Lb;
        ns_temp(I) = Lb(I);
        
        J = s>Ub;
        ns_temp(J) = Ub(J);
        
        s = ns_temp;
        
 %objective Function
 
        function z = fobj(u,prob)
            z = fitnessfuncsc208(u,prob);
%             if length(u) == 6
%                 z = fitnessfuncsc20(u,prob);
%             elseif length(u) == 15
%                 z = fitnessfuncsc205(u,prob);
%             end
