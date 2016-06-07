function [newpop]=insert_in_pop (cpopulation,nparticle,position)
newpop=[cpopulation(1:position-1,:);nparticle;cpopulation(position:end-1,:)];
