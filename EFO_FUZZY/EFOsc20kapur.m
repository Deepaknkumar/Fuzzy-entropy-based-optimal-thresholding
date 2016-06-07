%**************************************************************************************************
%Reference:  Abedinpourshotorban, H., Shamsuddin, S. M., Beheshti, Z., & Jawawi, D. N. (2015). 
%            Electromagnetic field optimization: A physics-inspired metaheuristic optimization algorithm. 
%            Swarm and Evolutionary Computation..
%
% Note: developed by Abedinpourshotorban, H. (2015)
%**************************************************************************************************

function [besterr,bestvalues,fminval]= EFOsc20kapur(N_var,N_emp,Max_gen,minval,maxval,R_rate,Ps_rate,P_field,N_field,E)
phi = (1 + sqrt(5))/2;%golden ratio
%initializatin
em_pop = minval + (maxval - minval).*rand(N_emp, N_var);
em_pop = sort(round(em_pop),2);

% fit= cec14_func(em_pop(1:N_emp,1:N_var)');  %%% fitness func
for i=1:N_emp
    fit(i) = fitnessfuncsc20kapur(em_pop(i,:),E); 
end

em_pop = [em_pop(:,1:N_var) fit'];
em_pop = sortpop(em_pop,N_var+1);

%random vectors (this is to increase the calculation speed instead of
%determining the random values in each iteration we allocate them in the
%beginning before algorithm start

r_index1 = randi([1 (round(N_emp.* P_field))],[N_var Max_gen]); %random particles from positive field
r_index2 = randi([(round(N_emp.*(1-N_field))) N_emp],[N_var Max_gen]); %random particles from negative field
r_index3 = randi([(round(N_emp.* P_field)+1) (round(N_emp.*(1-N_field))-1)],[N_var Max_gen]);%random particles from neutral field
ps= rand(N_var,Max_gen);%= probability of selecting electromagnets of generated particle from the positive field
r_force = rand(1,Max_gen);%random force in each generation
rp = rand(1,Max_gen);%some random numbers for checking randomness probability in each generation
randomization = rand(1,Max_gen);%coefficient of randomization when generated electro magnet is out of boundary
RI=1;%index of the electromagnet (variable) which is going to be initialized by random number
generation=N_emp;
new_emp = zeros(1,N_var+1); %temporary array to store generated particle 

while (generation <= Max_gen)
    r = r_force(1,generation);
    
    for i=1:N_var
              
        if (ps(i,generation) > Ps_rate)
            new_emp(i) = em_pop(r_index3(i,generation), i) + phi * r * (em_pop(r_index1(i,generation), i) - em_pop(r_index3(i,generation), i)) + r * (em_pop(r_index3(i,generation), i) - em_pop(r_index2(i,generation), i));
        else
            new_emp(i) = em_pop (r_index1(i,generation), i);
        end
        
        %checking whether the generated number is inside boundary or not
        if ( new_emp(i) >= maxval || new_emp(i) <= minval )
            new_emp(i) = minval + (maxval - minval) .* randomization(1, generation);
        end
    end
    new_emp = sort(round(new_emp(1:N_var)),2);
    %replacement of one electromagnet of generated particle with a random number (only for
    %some generated particles) to bring diversity to the population
    if ( rp(1,generation) < R_rate)
        new_emp(RI) = minval + (maxval - minval) .* randomization(1, generation);
        RI=RI+1;
        
        if (RI > N_var)
            RI=1;
        end
    end
    new_emp = sort(round(new_emp(1:N_var)),2);
%     new_emp(N_var+1) = cec14_func(new_emp(1:N_var)',problemIndex);
        new_emp(N_var+1) = fitnessfuncsc20kapur(new_emp(1:N_var),E); 
    %updating the population if the fitness of the generated particle is better than worst fitness in
    %the population (because the population is sorted by fitness, the last particle is the worst)
    if ( new_emp(N_var+1) < em_pop(N_emp , N_var+1) )
        position=find(em_pop(:,N_var+1) > new_emp(N_var+1));
        em_pop=insert_in_pop(em_pop,new_emp,position(1));
    end
    fminval(generation) = em_pop(1,N_var+1);
    generation=generation+1;
end
besterr=em_pop(1,N_var+1);
bestvalues = em_pop(1,1:N_var)
end