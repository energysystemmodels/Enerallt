% This file gets the information of network capacity and supply-demand curve of some interconnected regions and...
% determines the system price, minimizes area prices, and maximizes power flow in Nordpool (power market optimization)
% Copied from f_np_7R_dd
%----------------------------------------------------------------------------
function [Area_p,Sup_tot, Exch_opt,Ind_last,Sup_last,Sys_p]=f_np_dd(Demand, Pow0, Cost0,Ntc,Coff,Num_cost,Region_sys,peak_price)

% _________________________________________________________________________________
% 1) Preprations
n=size(Demand,1);
nh=size(Demand,2);
Demand=fix(Demand);

decim=10;                           % Making costs bigger to avoid decimals (optional)
Supply=cell(n,nh);                 % initial accumulative supply curve
Sup_dis=cell(n,nh);                % Supply curve after cost segmentation and moderation
%-----------------------------------------------------------------------
% 2) Lifting negative prices up, and later bringing them back
Cost_min=zeros(n,nh);
for t=1:nh
    for i=1:n
        Cost_min(i,t)=min(Cost0{i,t});
    end
end
cost_min=min(min(Cost_min));       % The minimum amount of cost in the whole 24 h

% Moderation of supply curve with the cost range and cost segments
Cost1=cell(n,nh);
Pow1=cell(n,nh);
for i=1:n
          for t=1:nh
              Supply{i,t}=cumsum(Pow0{i,t});
             
              [Sup_dis{i,t},Cost1{i,t}]=suppcurve_dis( Supply{i,t}, Cost0{i,t}-cost_min , Coff(i),Num_cost(i));
              Pow1{i,t}=[0 Sup_dis{i,t}(1) diff(Sup_dis{i,t})]; 
              Cost1{i,t}=[0 fix(decim.*(Cost1{i,t}))];     % MUST BE RETURNED LATER (Make them 10 times bigger and get rid of deciamls) 
          end

end
%-------------------------------------------------------------------------------------------------
% 2) Calculation of system price
% Crossing supply-demand curves to determine the system price (minimum possible price for meeting the whole demand)
my_subindex = @(A,r,c) A(r,c);      % An anonymous function to index a matrix (value = my_subindex(Matrix,row,column)

Sys_p=zeros(nh,1);
for t=1:nh
    N1=[];
    N2=[];
    for i=Region_sys
        N1=[N1 (Pow1{i,t})];
        N2=[N2 (Cost1{i,t})];
    end
    [Cost_sys, I_ns]=sort(N2);
    Pow_sys=N1(I_ns);
    Sys_p(t,1)=Cost_sys(find(cumsum(Pow_sys)>=sum(Demand(Region_sys,t),1),1));   % Initial system price
end
%-----------------------------------------------------------------------------------------------
% 3) Network and market optimization 
% The results are the optimal flow exchanges (Exch_opt) and final
% production in each region (Sup_tot), and cost of the whole system in euro
[Sup_tot,Exch_opt0,sys_cost]=f_opt_mar_dd(Demand,Cost1,Pow1,Ntc);
Exch_opt=Exch_opt0;
clear f_opt_mar_dd
%---------------------------------------------------------------------------------------------------------
% 4) Formation of area prices
Sys_p=Sys_p./decim+cost_min;
Area_p=zeros(n,nh);    % Area prices
for t=1:nh
    for i=1:n 
        if Sup_tot(i,t)> sum(Pow1{i,t})
            warning(['Region number ' num2str(i) ' is not capable to supply the domestic power without further import!'])
            Area_p(i,t)=peak_price;    % An arbitrary value to show the cost of reserve capacity activated
        elseif n==Region_sys
        Area_p(i,t)=max((Cost1{i,t}(find(((Sup_dis{i,t}))>=Sup_tot(i,t)-1,1))),Sys_p(t));
        else 
        Area_p(i,t)=(Cost1{i,t}(find(((Sup_dis{i,t}))>=Sup_tot(i,t)-1,1)));  
        end
    end
end

% Applying the Nordic power market rules: 
% 1. The importing country must have higher or eqaul price compared to the exporting
% 2. The price in region B which imports power from region A, is lower than all other areas exporting to region A.
for t=1:nh
    for i=1:n
        for j=1:n
            if Exch_opt(i,j,t)>0.1
                Area_p(j,t)=max([Area_p(j,t),Area_p(i,t),max(Area_p(:,t).*(Exch_opt(:,i,t)>0))]);
            end
        end
    end
end
    
Area_p=Area_p./decim+cost_min;
%-------------------------------------------------------------------------------------
% 4.1) Positioning of the last power producing unit in the power matrix

Ind_last=zeros(n,nh);      % The position of the last domestic producing plant in power production mix of that country
Sup_last=zeros(n,nh);      % Production from the last unit
for t=1:nh
    for i=1:n                % PAY ATTENTION
    Ind_last(i,t)=find(Supply{i,t}>=Sup_tot(i,t),1);
      if Ind_last(i,t)==1
        Sup_last(i,t)=Sup_tot(i,t);
      else
        Sup_last(i,t)=Sup_tot(i,t)-my_subindex(Supply{i,t},1,(Ind_last(i,t)-1));
      end
    end
end
%---------------------------------------------------------------------------------------------------------
% 5) Checking energy balances, transmission capacities, etc
for t=1:nh
    for i=1:n
        if (Demand(i)-sum(Exch_opt(:,i,t),1)-sum(Exch_opt(i,:,t),2))>1 && (Demand(i)-sum(Exch_opt(:,i,t),1)-sum(Exch_opt(i,:,t),2))<-1
         error(['Energy balance in the region number ' num2str(i) ' is not correct!!'])
        end
    end
    if abs(sum(Sup_tot(:,t))-sum(Demand(:,t)))>1 
        error('Energy balance in the whole system is not correct!!')
    end
    for i=1:n
        for j=1:n
          if (Ntc(i,j,t)-Exch_opt(i,j,t)<-1)
          error(['Power exchange in capacities higher than the interconnector between regions ' num2str(i) ' and ' num2str(j) ' !!']) 
          end
          if Area_p(i,t)>Area_p(j,t)+0.5 && Exch_opt(i,j,t)>1
          warning(['Power transmission from a higher area price to a lower area price, from region ' num2str(i) ' to region ' num2str(j) ': ' num2str(Exch_opt(i,j,t)) ' MW !!'])
          end
        end
    end
    for i=1:n
        for j=1+i:n
            if Exch_opt(i,j,t)>1 && Exch_opt(j,i,t)>1
                warning([num2str(Exch_opt(i,j,t)) ' and ' num2str(Exch_opt(j,i,t))])
               %error(['Power transmission in both directions in the line between regions ' num2str(i) ' and ' num2str(j)' ' !!'])
            end
        end
    end
end