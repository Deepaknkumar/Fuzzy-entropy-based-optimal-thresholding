%% MODE
% Multi-objective Evolutionary Algorithm (MOEA) based on Differential
% Evolution (DE).
% It implements a greedy selection based on pure dominance.
% DE algorithm has been introduced in:
%
% Storn, R., Price, K., 1997. Differential evolution: A simple and 
% efficient heuristic for global optimization over continuous spaces. 
% Journal of Global Optimization 11, 341 – 359.
%%

function [threshvalues,minval,mud,mum,mub,bestnest] =DEsc20(num,number_levels,E)

NOBJ = 1;                          % Number of objectives
MODEDat.NRES = 0;                          % Number of constraints
Nvar   = number_levels*3;       % Numero of decision variables
nub = Nvar;
MODEDat.mop = str2func('fitnessfunc');    % Cost function

MODEDat.FieldD =[zeros(nub,1) 255.*ones(nub,1)]; % Initialization bounds
MODEDat.Initial=[zeros(nub,1) 255.*ones(nub,1)]; % Optimization bounds

%MODEDat.XPOP = n*NOBJ;             % Population size
% MODEDat.Esc = 0.5;                         % Scaling factor
% MODEDat.Pm= 0.2;                           % Croosover Probability

%MODEDat.MAXGEN =10000;                     % Generation bound
MODEDat.MAXFUNEVALS = 150*Nvar...  % Function evaluations bound
    *NOBJ;   

MODEDat.CounterGEN=0;
MODEDat.CounterFES=0;
%% Reading parameters from MODEDat
Generaciones  = 700; %MODEDat.MAXGEN;    % Maximum number of generations.
Xpop          =   num*NOBJ; %MODEDat.XPOP;      % Population size.
%Nvar;          %= Nvar;      % Number of decision variables.
Nobj          = 1; % MODEDat.NOBJ;      % Number of objectives.
Bounds        = MODEDat.FieldD;    % Optimization bounds.
Initial       = MODEDat.Initial;   % Initialization bounds.
ScalingFactor = 0.5;  %MODEDat.Esc;       % Scaling fator in DE algorithm.
CrossOverP    = 0.2;  %MODEDat.Pm;        % Crossover probability in DE algorithm.
mop           = MODEDat.mop;       % Cost function.
MODEDat.InitialPop = [];


%% Initial random population
Parent = zeros(Xpop,Nvar);  % Parent population.
Mutant = zeros(Xpop,Nvar);  % Mutant population.
Child  = zeros(Xpop,Nvar);  % Child population.
FES    = 0;                 % Function Evaluation.

for xpop=1:Xpop
    for nvar=1:Nvar
        Parent(xpop,nvar) = Initial(nvar,1)+(Initial(nvar,2)...
                            - Initial(nvar,1))*rand();
    end
end

Parent = sort(round(Parent),2);

if size(MODEDat.InitialPop,1)>=1
    Parent(1:size(MODEDat.InitialPop,1),:)=MODEDat.InitialPop;
end

sof = size(Parent);
for i=1:sof(1)
    JxParent(i,1) = fitnessfuncsc20(Parent(i,:),E);
end
FES = FES+Xpop;   

%% Evolution process

for n=1:Generaciones 
    
    for xpop=1:Xpop
        rev=randperm(Xpop);
        
        %% Mutant vector calculation
        Mutant(xpop,:)= Parent(rev(1,1),:)+ScalingFactor*...
                       (Parent(rev(1,2),:)-Parent(rev(1,3),:));
        Mutant = sort(round(Mutant),2);
        for nvar=1:Nvar %Bounds are always verified
            if Mutant(xpop,nvar)<Bounds(nvar,1)
                Mutant(xpop,nvar) = Bounds(nvar,1);
            elseif Mutant(xpop,nvar)>Bounds(nvar,2)
                Mutant(xpop,nvar)=Bounds(nvar,1);
            end
        end
        
        %% Crossover operator
        for nvar=1:Nvar
            if rand() > CrossOverP
                Child(xpop,nvar) = Parent(xpop,nvar);
            else
                Child(xpop,nvar) = Mutant(xpop,nvar);
            end
        end

    end
    Child = sort(round(Child),2);
    sofc = size(Child);
for i=1:sofc(1)
    JxChild(i,1) = fitnessfuncsc20(Child(i,:),E);
end

    FES=FES+Xpop;

    %% Selection
    for xpop=1:Xpop
        if JxChild(xpop,:) <= JxParent(xpop,:) 
            Parent(xpop,:) = Child(xpop,:);
            JxParent(xpop,:) = JxChild(xpop,:);
        end
    end
    
	PFront=JxParent;
	PSet=Parent;

    OUT.Xpop           = Parent;   % Population
    OUT.Jpop           = JxParent; % Poopulation's Objective Vector
    OUT.PSet           = PSet;     % Pareto Set
    OUT.PFront         = PFront;   % Pareto Front
    OUT.Param          = MODEDat;  % MODE Parameters
    MODEDat.CounterGEN = n;
    MODEDat.CounterFES = FES;
    
%     [OUT MODEDat]=PrinterDisplay(OUT,MODEDat); % To print results on screen
%     
%     if FES>MODEDat.MAXFUNEVALS || n>MODEDat.MAXGEN
%         disp('Termination criteria reached.')
%         break;
[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet);

%  bnvalues(n,:) = OUT.PSet;
end


OUT.Xpop=PSet;
OUT.Jpop=PFront;
[OUT.PFront, OUT.PSet]=DominanceFilter(PFront,PSet); %A Dominance Filter
u = OUT.PSet;
minval = OUT.PFront;

bestnest = u
threshvalues = [];
if u(2) >= (u(1)+u(3))/2
    threshvalues(1) = u(1) + sqrt((u(3)-u(1))*(u(2)-u(1))/2);
elseif u(2) < (u(1)+ u(3))/2
    threshvalues(1) = u(3) - sqrt((u(3)-u(1))*(u(3)-u(2))/2);
end

if u(5) >= (u(4)+u(6))/2
    threshvalues(2) = u(4) + sqrt((u(6)-u(4))*(u(5)-u(4))/2);
elseif u(5) < (u(4)+u(6))/2
    threshvalues(2) = u(6) - sqrt((u(6)-u(4))*(u(6)-u(5))/2);
end
     
            mud = zeros(1,256);
            mum = zeros(1,256);
            mub = zeros(1,256);
            for k=1:256
                if k<=u(1)
                    mud(k) = 1;
                    mum(k) = 0;
                    mub(k) = 0;
                end
                if k<=u(2) && k>u(1)
                    mud(k) = 1 - (((k-u(1))^2)/((u(3)-u(1))*(u(2)-u(1))));
                    mum(k) = (((k-u(1))^2)/((u(3)-u(1))*(u(2)-u(1))));
                    mub(k) = 0;
                end
                if k>u(2) && k<=u(3)
                    mud(k) = (((k-u(3))^2)/((u(3)-u(1))*(u(3)-u(2))));
                    mum(k) = 1 - (((k-u(3))^2)/((u(3)-u(1))*(u(3)-u(2))));
                    mub(k) = 0;
                end
                if k>u(3) && k<=u(4)
                    mud(k) = 0;
                    mum(k) = 1;
                    mub(k) = 0;
                end
                if k>u(4) && k<=u(5)
                    mud(k) = 0;
                    mum(k) = 1 - (((k-u(4))^2)/((u(6)-u(4))*(u(5)-u(4))));
                    mub(k) = (((k-u(4))^2)/((u(6)-u(4))*(u(5)-u(4))));
                end
                if k>u(5) && k<=u(6)
                    mud(k) = 0;
                    mum(k) = (((k-u(6))^2)/((u(6)-u(4))*(u(6)-u(5))));
                    mub(k) = 1 - (((k-u(6))^2)/((u(6)-u(4))*(u(6)-u(5))));
                end
                if k> u(6)
                    mud(k) = 0;
                    mum(k) = 0;
                    mub(k) = 1;
                end
            end


% if strcmp(MODEDat.SaveResults,'yes')
%     save(['OUT_' datestr(now,30)],'OUT'); %Results are saved
% end
% 
% disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
% disp('Red  asterisks : Set Calculated.')
% disp('Black diamonds : Filtered Set.')
% if strcmp(MODEDat.SaveResults,'yes')
%     disp(['Check out OUT_' datestr(now,30) ...
%           ' variable on folder for results.'])
% end
% disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
% 
% F=OUT.PFront;
% for xpop=1:size(F,1)
%     if Nobj==1
%         figure(123); hold on;
%         plot(MODEDat.CounterGEN,log(min(F(:,1))),'dk', ...
%             'MarkerFaceColor','k'); grid on; hold on;
%     elseif Nobj==2
%         figure(123); hold on;
%         plot(F(xpop,1),F(xpop,2),'dk','MarkerFaceColor','k');...
%             grid on; hold on;
%     elseif Nobj==3
%         figure(123); hold on;
%         plot3(F(xpop,1),F(xpop,2),F(xpop,3),'dk','MarkerFaceColor','k');...
%             grid on; hold on;
%     end
% end

%% Print and Display information
% Modify at your convenience
%
% function [OUT Dat]=PrinterDisplay(OUT,Dat)
% 
% disp('------------------------------------------------')
% disp(['Generation: ' num2str(Dat.CounterGEN)]);
% disp(['FEs: ' num2str(Dat.CounterFES)]);
% disp(['Pareto Front Size: ' mat2str(size(OUT.PFront,1))]);
% disp('------------------------------------------------')
% 
% if mod(Dat.CounterGEN,1)==0
%     if Dat.NOBJ==3
%         figure(123);
%         plot3(OUT.PFront(:,1),OUT.PFront(:,2),OUT.PFront(:,3),'*r'); 
%         grid on;
%     elseif Dat.NOBJ==2
%         figure(123);
%         plot(OUT.PFront(:,1),OUT.PFront(:,2),'*r'); grid on;
%     elseif Dat.NOBJ==1
%         figure(123);
%         plot(Dat.CounterGEN,log(min(OUT.PFront(:,1))),'*r'); ...
%             grid on; hold on;
%     end
% end

%% Dominance Filter
%
% A filter based on dominance criteria
%
function [Frente Conjunto]=DominanceFilter(F,C)

Xpop=size(F,1);
Nobj=size(F,2);
Nvar=size(C,2);
Frente=zeros(Xpop,Nobj);
Conjunto=zeros(Xpop,Nvar);
k=0;

for xpop=1:Xpop
    Dominado=0;
    
    for compara=1:Xpop
        if F(xpop,:)==F(compara,:)
            if xpop > compara
                Dominado=1;
                break;
            end
        else
            if F(xpop,:)>=F(compara,:)
                Dominado=1;
                break;
            end
        end
    end
    
    if Dominado==0
        k=k+1;
        Frente(k,:)=F(xpop,:);
        Conjunto(k,:)=C(xpop,:);
    end
end
Frente=Frente(1:k,:);
Conjunto=Conjunto(1:k,:);

%% Release and bug report:
%
% November 2012: Initial release