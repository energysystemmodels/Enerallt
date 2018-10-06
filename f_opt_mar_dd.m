% Optimization Module
% This function receives the supply-cost curve of production modes in several nodes, 
% network capacities and costs, and demand at each node. This function optmizes the 
% flow through the network to minimize the total operation cost of the whole system.
%__________________________________________________________________________
% 1) Input arguments
% Choose only one of the three methods below

% 1.1) Using this module as a function
%function [Sup_tot, Exch_opt, sys_cost]=f_opt_mar_dd(Demand,Cost,Pow,NTC)

%--------------------------------------------------------------------------
% 1.2) Reading input data from this file
tic
h=1;
Demand=[450 800 600 350]';      % power demand in each node
Demand=repmat(Demand,1,h);
n=size(Demand,1);

P1=[200 300 300];
P2=[140 300 100 200];
P3=[240 100 250];
P4=[150 100 200 100 200];

C1=[30 41 82];
C2=[44 45 46 80];
C3=[30 60 80];
C4=[34 44 50 10 10];


Coff=[0.2 0.2 .2 0.2];
% from (row) and to (coloumn)
Ntc=[0   60  90   50;
     60  0   50   0;
     90  50  0    100;
     0   0   100  0];
NTC=repmat(Ntc,1,1,h);


Pow=cell(n,h);
Cost=cell(n,h);

for t=1:h
Pow{1,t}=P1;Pow{2,t}=P2;Pow{3,t}=P3;Pow{4,t}=P4;
Cost{1,t}=C1;Cost{2,t}=C2;Cost{3,t}=C3;Cost{4,t}=C4;
end

%--------------------------------------------------------------------------
% 1.3) reading input data from workspace (from f_nd_7R_dd)
% vis=0;
% Demand=evalin( 'base', 'Demand');
% Pow=evalin( 'base', 'Pow1');
% Cost=evalin( 'base', 'Cost1');
% NTC=evalin( 'base', 'Ntc');
% %cost_min=evalin( 'base', 'cost_min');
% Coff=evalin( 'base', 'Coff');

%____________________________________________________________________________
% 2) Preparation

n=size(Demand,1);
h=size(Demand,2);

% Making a matrix with the same length of the costs
Cost_min=zeros(n,h);
sumpow=zeros(n,h);
for t=1:h
    for i=1:n
        Cost_min(i,t)=min(Cost{i,t});         % Finding minimum cost
        sumpow(i,t)=sum(Pow{i,t});
    end
    if sum(sumpow(:,t))<sum(Demand(:,t))      % Check if there is enough supply in the system
            error(['There is no enough power supply to meet demand at time step h = ' num2str(t) ' !'])
    end
end
cost_min=min(min(Cost_min)); 
%__________________________________________________________________________
% 3) Modeling the problem (setting up the optimization)

% Variables: electricity supply from each production mode in each region + 
% electricity trade from each node to the interconnected nodes

% TODO: At the moment, the trade possibilites between all nodes are considered 
% as variables. However, it should only consider those trade options that have 
% network capacity (using "find" command). This should reduce the number of variables.


% 3.1) Objective function
f=[];
for t=1:h        % Number of timesteps in each optimization run (e.g., h=24 for hourly optimization over one day)
    for i=1:n    % Number of nodes
    f=[f repmat((Cost{i,t}),1,1+n)];              % Cost function, repeated n+1 times (one time for power supply internally, and n times to all nodes) 
    end
end
f=f-cost_min;                  % The minimum amount of cost in the whole period is used for adjusting all the costs positive

%------------------------------------
% 3.2) Electricity transmission cost
t_cost = 0.0001;              % It can be based on each line length and capacity
nvar=0;                       % Counter of decision variables (supply modes in each node multiplied by (n+1))
for t=1:h                     % Time steps in each run
    for i=1:n
      for j=1:n
          f(nvar+j*length(Cost{i,t})+1:nvar+(j+1)*length(Cost{i,t})) = f(nvar+j*length(Cost{i,t})+1:nvar+(j+1)*length(Cost{i,t})) + t_cost;
      end
    nvar=length(Cost{i})*(n+1)+nvar;
    end
end
%---------------------------------------
% 3.3) Bounds
NTC_rep=cell(n,h);           % Maximum network capacity to be reshaped to match objective function
for t=1:h
    for i=1:n
        NTC_rep{i,t}=reshape(permute(repmat(NTC(i,:,t)',1,length(Cost{i,t})),[2 1]),1,n*length(Cost{i,t}));
    end
end

% 3.3.1) Upper bounds
UB_1=cell(n,h);
UB=[];
for t=1:h
    for i=1:n
    UB_1{i,t}=[Pow{i,t} NTC_rep{i,t}];           
    UB=[UB UB_1{i,t}];
    end
end

% 3.3.2) Lower bounds
LB=zeros(size(f));

%-----------------------------------------
% 3.4) Constraints
% 3.4.1) Inequality constraints (A.X =< b)
A=[];
b=[];

% The sum of electricity export in each node to other nodes shouldn't exceed the max network capacity (n^2*h equations)
A1=spalloc(n^2*h,length(f),(n+1)*length(Cost{1,1})*h*n); 
b1=zeros(n^2*h,1);              % Matrix of consonents  

hrow=0;                         % Counter of row blocks (each block for one time step h)
hcol=0;                         % Counter of column blocks (each block for one time step h) 
for t=1:h
    nnode=0;         % Counter of node blocks (shifting rows)
    nvar=0;          % Counter of decision variables (supply modes in each node multiplied by (n+1)) 
    ncost=0;         % Counter of supply modes in each node in each time step
    for i=1:n
        for j=1:n
            A1(hrow+nnode+j, hcol+nvar+j*length(Cost{i,t})+1 : hcol+nvar+(j+1)*length(Cost{i,t}))=1;     
        end
         nvar=length(Cost{i,t})*(n+1)+nvar;
         nnode=n+nnode;
         ncost=length(Cost{i,t})+ncost;
    end
    b1(hrow+1:hrow+n*n)=reshape(permute(NTC(:,:,t),[2 1 3]),n^2,1);
    hrow=hrow+n*n;
    hcol=hcol+nvar;  
end

ncost_tot = ncost;
% The sum of power production and export from each supply mode shouldn't exceed the max power capacity of that mode
A2= spalloc(ncost_tot*h,length(f),n*length(Cost{1,1})*h*n);   % Number of equations = supply modes in all nodes multiplied by time steps)
b2=[];
hrow=0;
hcol=0;
for t=1:h
    ncost=0;
    nvar=0;
    for i=1:n
            A2(hrow+1+ncost:hrow+ncost+length(Cost{i,t}),hcol+nvar+1:hcol+nvar+(n+1)*length(Cost{i,t}))=repmat(eye(length(Cost{i,t})),1,n+1);
            ncost=ncost+length(Cost{i,t});
            nvar=length(Cost{i,t})*(n+1)+nvar;
            b2=[b2 Pow{i,t}];
    end
    hrow=hrow+ncost_tot;
    hcol=hcol+nvar;
end

b2=b2'; 
A=[A1;A2];
b=[b1;b2];

% 3.4.2) Equality constraints (Aeq.X = beq)
% Energy balance for each node (n equations)
Aeq1=spalloc(n*h,length(f),n*length(Cost{1,1})*h*n);

hrow=0;
hcol=0;
for t=1:h
    nvar=0;
    for i=1:n
        Aeq1(hrow+i,hcol+nvar+1:hcol+nvar+length(Cost{i,t}))=1;          % Production of own plants
        Aeq1(hrow+i,hcol+nvar+length(Cost{i,t})+1:hcol+nvar+(n+1)*length(Cost{i,t}))=-1;     % Import from other regions

        %Aeq1(hrow+i,hcol+z+i*length(Cost{i,t})+1:hcol+z+(i+1)*length(Cost{i,t}))=1;       % Just for prepration for the following loop
        %Aeq2(hrow+z+1:hrow+z+length(Cost{i,t}))=1;                             % Production of plants equal to demand
        nvar2=0;
        for j=1:n 
            Aeq1(hrow+i,hcol+nvar2+i*length(Cost{j,t})+1:hcol+nvar2+(i+1)*length(Cost{j,t}))=Aeq1(hrow+i,hcol+nvar2+i*length(Cost{j,t})+1:hcol+nvar2+(i+1)*length(Cost{j,t}))+1; % Export to other regions
            nvar2=length(Cost{j,t})*(n+1)+nvar2;
        end
        nvar=length(Cost{i,t})*(n+1)+nvar;
    end
    hrow=hrow+n;
    hcol=hcol+nvar;
end

beq1=reshape(Demand(:,1:h),n*h,1);                  % A column vector in which the number of rows denote the number of equality constraints

%__________________________________________________________________________
% 4) Solving the model and postprocessing
% 4.1) Optimization
options = optimoptions(@linprog,'Display','off','MaxIter',20000,'Algorithm','dual-simplex','TolFun',1e-6,'TolCon',1e-6);
[Opt ,sys_cost0]=linprog(f,A1,b1,Aeq1,beq1,LB,UB,[],options);

if isempty(Opt) == 1
    warning('No optimal solution found, zero exchange!!')
    Opt = zeros(size(f))';
end
sys_cost=sum(Opt'.*(f+cost_min));   % The corrected system cost: Deducting the value of cost_min from cost items, which was added before to handle negative values
Exch_opt=zeros(size(NTC));          % Optimal trade flow
Pow_final=cell(n,h);                % Final supply by each production mode (accepted bids)
Sup_tot=zeros(n,h);                 % Total production in each node (incl. export to other nodes)

nvar=0;
sumpp=0;
for t=1:h
    for i=1:n
        Pow_final{i,t}=Opt(nvar+1:nvar+length(Cost{i,t}));
        for j=1:n
            Exch_opt(i,j,t)=sum(Opt(nvar+j*length(Cost{i,t})+1:nvar+(j+1)*length(Cost{i,t})));
        end
        nvar=length(Cost{i,t})*(n+1)+nvar;
        sumpp=sumpp+sum(Pow_final{i,t});
    end
    Sup_tot(:,t)=Demand(:,t)-sum(Exch_opt(:,:,t),1)'+sum(Exch_opt(:,:,t),2);
end

en_balance=fix(sum(sum(Demand))-sum(sum((Sup_tot))));               % Test: energy balance
Exch_opt=fix(Exch_opt);
Sup_tot=fix(Sup_tot);
%--------------------------------------------------------------------------
% 4.2) Postprocessing
% 4.2.1) Area prices
Supp_curve=cell(n,1);
for i=1:n
    [Cost{i}, I]=sort(Cost{i});
    Pow{i}=Pow{i}(I);
    Supp_curve{i}=suppcurve( cumsum(Pow{i}), Cost{i}, Coff(i));
end

Area_p=zeros(n,1);
for i=1:n
    if fix(Sup_tot(i))==0
        warning(['No need for power production in region '  num2str(i)  ' !*'])
        Area_p(i)=[];
    else
    Area_p(i)=Supp_curve{i}(Sup_tot(i));
    end
end

toc
                                   
%--------------------------------------------------------------------------
% 4.2.2) Area prices before the trade
my_subindex = @(A,r,c) A(r,c);      % An anonymous function to index a matrix
                                    % value = my_subindex(Matrix,row,column); 

% 4.2.3) System price 
P=[];
C=[];
for i=1:n
P=[P (Pow{i})];
C=[C (Cost{i})];
end
[Cost_sys, I_ns]=sort(C);
Pow_sys=P(I_ns);
                                                                      
Supp_cur_sys=suppcurve( cumsum(Pow_sys), Cost_sys, mean(Coff));
sys_p=my_subindex(Supp_cur_sys,1,(find(1:max(cumsum(Pow_sys))>=sum(Demand),1)));  % Initial system price
%__________________________________________________________________________
% 5) Visualization
figure
Area_p0=zeros(n,1);
for i=1:n
subplot(round(n/2+1),2,i)
plot(1:max(cumsum(Pow{i})),Supp_curve{i});hold on;
% Before trade
if max(cumsum(Pow{i}))<Demand(i)
    warning(['Region number ' num2str(i) ' is not able to supply its demand without imports!'])
    Area_p0(i)=250;
else
    Area_p0(i)=my_subindex(Supp_curve{i},1,Demand(i));
end
X=[0 Demand(i) ;Demand(i) Demand(i)];
Y=[Area_p0(i) Area_p0(i) ; Area_p0(i) 0 ];
line(X,Y,'Color','r','LineStyle','--');hold on;
% After trade
X=[0  Sup_tot(i);Sup_tot(i) Sup_tot(i)];
Y=[Supp_curve{i}(Sup_tot(i)) Supp_curve{i}(Sup_tot(i)); Supp_curve{i}(Sup_tot(i)) 0 ];
line(X,Y,'Color','black','LineStyle','-.');hold on;
xlabel('Quantity (MW)');ylabel('Price (€/MWh)');title(['Node ' num2str(i)]);
xlim([0 1.1*max(cumsum(Pow{i}))])
ylim([0 1.1*max(Supp_curve{i})])
text(0.2*Demand(i) ,Area_p(i)+4,[num2str(Area_p(i) ,'%.1f') '€'],'FontSize',8);
end

subplot(round(n/2+1),2,i+1)
plot(1:max(cumsum(Pow_sys)),Supp_cur_sys);hold on;
X=[0 sum(Demand); sum(Demand) sum(Demand)];
Y=[sys_p sys_p; sys_p 0];
line(X,Y,'Color','r','LineStyle','--');hold on;
xlabel('Quantity (MW)');ylabel('Price (€/MWh)');title('Nordic system');
xlim([0 1.1*max(cumsum(Pow_sys))])
text(0.2*sum(Demand) ,sys_p+4,[num2str(sys_p ,'%.2f') '€'],'FontSize',8); 

