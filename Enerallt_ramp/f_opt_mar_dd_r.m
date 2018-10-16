% Optimization Module
% This function receives the supply-cost curve of production modes in several nodes, 
% network capacities and costs, and demand at each node. This function optmizes the 
% flow through the network to minimize the total operation cost of the whole system.
%__________________________________________________________________________
% 1) Input arguments
function [Sup_tot, Exch_opt, sys_cost]=f_opt_mar_dd_r(Demand,Pow,Cost,Ramp,NTC)
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
        nvar=length(Cost{i,t})*(n+1)+nvar;
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

    prod_mode_tot = ncost;
    %--------------------------------------------------------------------------------------
    % 3.4.1.2) Ramping constraints
    % Upward ramping: x2 =< x1 + ramp*x1 ==> x2 - (1+ramp)*x1 =< 0
    A2= spalloc(prod_mode_tot*(h-1),length(f),n*length(Cost{1,1})*h*n);   % Number of equations = production modes in all nodes multiplied by time steps minus one (no constraint for the first hour)
    b2=zeros(prod_mode_tot*(h-1),1);
    hrow=0;
    hcol=0;
    for t=1:h-1
        ncost=0;
        nvar=0;
        for i=1:n
            A2(hrow+ncost+1:hrow+ncost+length(Cost{i,t}),hcol+nvar+1:hcol+nvar+length(Cost{i,t}))...
                =-eye(length(Cost{i,t}));   % Fisrt hour
            for j=1:length(Cost{i,t+1})     % Consecutive hour
                A2(hrow+ncost+j,hcol+nvar+prod_mode_tot*(n+1)+j)=1/(1+Ramp{i,t+1}(j));   
                b2(hrow+ncost+j,1)=0.33*Pow{i,t+1}(j)+300;
            end 
                ncost=ncost+length(Cost{i,t});
                nvar=length(Cost{i,t})*(n+1)+nvar;
        end
        hrow=hrow+prod_mode_tot;
        hcol=hcol+nvar;
    end 

    % Downward ramping: x2 >= x1 - ramp*x1 ==> (1-ramp)*x1 - x2 =< 0
    % A3= spalloc(prod_mode_tot*(h-1),length(f),n*length(Cost{1,1})*h*n);   % Number of equations = production modes in all nodes multiplied by time steps minus one (no constraint for the first hour)
    % b3=zeros(prod_mode_tot*(h-1),1);
    % hrow=0;
    % hcol=0;
    % for t=1:h-1
    %     ncost=0;
    %     nvar=0;
    %     for i=1:n
    %         A3(hrow+ncost+1:hrow+ncost+length(Cost{i,t}),hcol+nvar+1:hcol+nvar+length(Cost{i,t}))...
    %             =eye(length(Cost{i,t}));   % First hour
    %         for j=1:length(Cost{i,t+1})
    %             A3(hrow+ncost+j,hcol+nvar+prod_mode_tot*(n+1)+j)=-1/(1-Ramp{i,t+1}(j));   % Consecutive hour
    %             b3(hrow+ncost+j,1)=0.33*Pow{i,t+1}(j)+300;      # SHOULD BE THOUGHT MORE
    %         end    
    %             ncost=ncost+length(Cost{i,t});
    %             nvar=length(Cost{i,t})*(n+1)+nvar;
    %     end
    %     hrow=hrow+prod_mode_tot;
    %     hcol=hcol+nvar;
    % end
    %---------------------------------------------------------------------------
    % 3.4.2) Equality constraints (Aeq.X = beq)
    % Energy balance for each node (n equations)
    A4=spalloc(n*h,length(f),n*length(Cost{1,1})*h*n);

    hrow=0;
    hcol=0;
    for t=1:h
        nvar=0;
        for i=1:n
            A4(hrow+i,hcol+nvar+1:hcol+nvar+length(Cost{i,t}))=1;          % Production of own plants
            A4(hrow+i,hcol+nvar+length(Cost{i,t})+1:hcol+nvar+(n+1)*length(Cost{i,t}))=-1;     % Import from other regions

            %Aeq1(hrow+i,hcol+z+i*length(Cost{i,t})+1:hcol+z+(i+1)*length(Cost{i,t}))=1;       % Just for prepration for the following loop
            %Aeq2(hrow+z+1:hrow+z+length(Cost{i,t}))=1;                             % Production of plants equal to demand
            nvar2=0;
            for j=1:n 
                A4(hrow+i,hcol+nvar2+i*length(Cost{j,t})+1:hcol+nvar2+(i+1)*length(Cost{j,t}))=...
                    A4(hrow+i,hcol+nvar2+i*length(Cost{j,t})+1:hcol+nvar2+(i+1)*length(Cost{j,t}))+1; % Export to other regions
                nvar2=length(Cost{j,t})*(n+1)+nvar2;
            end
            nvar=length(Cost{i,t})*(n+1)+nvar;
        end
        hrow=hrow+n;
        hcol=hcol+nvar;
    end

    b4=reshape(Demand(:,1:h),n*h,1);                  % A column vector in which the number of rows denote the number of equality constraints
    % Collecting all inaquality constraints
    A=[A1;A2;-A4];
    b=[b1;b2;-b4];
    %__________________________________________________________________________
    % 4) Solving the model and postprocessing
    % 4.1) Optimization
    options = optimoptions(@linprog,'Display','off','MaxIter',20000,'Algorithm','dual-simplex','TolFun',1e-4,'TolCon',1e-6);
    [Opt ,sys_cost0]=linprog(f,A,b,[],[],LB,UB,[],options);

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
    %en_balance=fix(sum(sum(Demand))-sum(sum((Sup_tot))))               % Test: energy balance
