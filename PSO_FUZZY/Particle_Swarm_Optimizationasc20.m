function [ gBest_availability,optimised_parameters,fminval] = Particle_Swarm_Optimizationasc20(Bird_in_swarm, Number_of_quality_in_Bird, MinMaxRange, Food_availability, availability_type, velocity_clamping_factor, cognitive_constant, social_constant, Min_Inertia_weight, Max_Inertia_weight, max_iteration,E)


if nargin < 11
    error('Missing input parameter(s)!')
end

%{
	universalize availability type  
%}
availability_type=lower(availability_type(1:3));

[row,col]=size(MinMaxRange);
if row~=Number_of_quality_in_Bird || col~=2
    error('Not a proper MinMaxRange Matrix')
end
for i=1:Number_of_quality_in_Bird
    if MinMaxRange(i,1)>=MinMaxRange(i,2)
        error('Minimum value greater than Maximum value!!!')
    end
end

%{
	 counter to display % of completion
%}
N=Bird_in_swarm*max_iteration;
q=0;

%{
	 distinguishing min and max range
%}
bird_min_range=MinMaxRange(:,1);
bird_max_range=MinMaxRange(:,2);

%{
	 
%}
format long;
for i=1:Number_of_quality_in_Bird
    bird(:,i)=bird_min_range(i)+(bird_max_range(i)-bird_min_range(i))*rand(Bird_in_swarm,1);
end
bird = sort(round(bird),2);

Vmax=bird_max_range*velocity_clamping_factor;
Vmin=-Vmax;

for i=1:Number_of_quality_in_Bird
    Velocity(:,i)=Vmin(i)+(Vmax(i)-Vmin(i))*rand(Bird_in_swarm,1);
end

for itr=1:max_iteration
%     fprintf('Completed  %d  %% ...', uint8(q*100/N ))
    
    for p=1:Bird_in_swarm
        parameter=bird(p,:,itr);
        availability(p,itr)=fitnessfuncsc208(parameter,E);
        
        switch availability_type
            case 'min'
                format long;
                [pBest_availability,index]=min(availability(p,:));
                pBest=bird(p,:,index);
                
                if(p==1 && itr==1)
                    gBest=pBest;
                    gBest_availability=pBest_availability;
                elseif availability(p,itr)<gBest_availability
                    gBest_availability=availability(p,itr);
                    gBest=bird(p,:,itr);
                end
                
            case 'max'
                format long;
                [pBest_availability,index]=max(availability(p,:));
                pBest=bird(p,:,index);
                
                if(p==1 && itr==1)
                    gBest=pBest;
                    gBest_availability=pBest_availability;
                elseif availability(p,itr)>gBest_availability
                    gBest_availability=availability(p,itr);
                    gBest=bird(p,:,itr);
                end
                
            otherwise
                error('availability_type mismatch')
        end
        
        w(itr)=((max_iteration - itr)*(Max_Inertia_weight - Min_Inertia_weight))/(max_iteration-1) + Min_Inertia_weight;
        Velocity(p,:,(itr+1))=w(itr)*Velocity(p,:,itr) + social_constant*rand(1,Number_of_quality_in_Bird).*(gBest-bird(p,:,itr)) + cognitive_constant*rand(1,Number_of_quality_in_Bird).*(pBest-bird(p,:,itr));
        Velocity(p,:,(itr+1))=MinMaxCheck(Vmin, Vmax, Velocity(p,:,(itr+1)));
        
        bird(p,:,(itr+1))= bird(p,:,itr) + Velocity(p,:,(itr+1));
        bird(p,:,(itr+1))=MinMaxCheck(bird_min_range, bird_max_range, bird(p,:,(itr+1)));
        q=q+1;
    end
        bird = sort(round(bird),2);
        fminval(itr, :) = gBest_availability;
end
gBest_availability;
optimised_parameters=gBest

threshvalues = [];
u = optimised_parameters;