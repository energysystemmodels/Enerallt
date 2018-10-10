% This model optimizes a power market and network between a group of interconnected energy systems
% This script requires the following functions to run
% 1. eP_function  % Reads the input data and simulates the energy systems
% 2. suppcurve    % Modifies supply cruve (and adjusts some increments if needed)
% 3. hydroplan    % Offering hydro based on water value concept
% 4. f_opt_sys    % Calculation of system price
% 5. f_opt_mar    % Market optimization module
% 6. f_np_6RD     % A facilitator for establishing the area prices after market optimization (and visualization if neeed)

%_______________________________________________________________________
% 1) Input data and desirable outputs
clear all
close all
save_name='Ref_7R';         % Name of the file for saving the results
Change=2;
h=24*4+10;                          % An hour to see the sample table results
%---------------------------------------------------------------------------------
% 1.1) Input data
td_1=1;                % The beginning day of the analysis
td_n=60;               % The last day of the analysis
n=6;                   % Number of regions to be simulated
Region_sys=1:5;        % The regions that are in one system price calculations
Region_wv=[1 2 5];     % Regions that have hydropower water value calculations
peak_price=1000;        % The electricity price when there is no enough suply in the market (to indicte those hours)
yr=2014;               % The year of the analysis
% t_d_fi=+1;           % Time difference between Finland and other countries (central EU +...)
%------------------------------------------------------------------------
% 2) Prepration and initialization 
% 2.1) Timeframe

nd=td_n-td_1+1;        % Number of days (if you wish for hours less than one day, you may consider changing this)

th_1=24*(td_1-1)+1;    
th_n=24*td_n;
if th_1>th_n
    error('Number of hours is smaller than at least one whole day!')
end
nh=th_n-th_1+1;        % Number of hours in the analysis
%----------------------------------------------------------------------------

[Ntc,Market_exch,El_dem_tot,El_dem_res_tot,Eff_el,Emm_fac,Eff_el_chp1,Eff_el_chp2,Eff_chp_boil,Eff_hob,....
    El_pro_vre1,El_pro_vre2,El_pro_vre3,El_pro_vre4,El_pro_nuc,El_chp_ind_hh,El_pro_fix_hh,...
    El_cap_hyd,El_eff_hyd,Hyd_res_dd,Hyd_infl_dd,St_hyd,Cost_hyd_corr,Wv_seg,Cost_wv,Cost_hyd_coff,St_begin,St_min,Ind_hyd,...
    El_pro_chp1,El_pro_chp,Pth_chp,Dh_pro_chp,Dh_pro_boil,Dh_pro_hob,F_chp1_r,F_chp2_r,F_boil_r,F_hob_r,El_chp1_mar,El_chp2_mar,El_ppchp2_mar,El_ppchp_mar,...
    Cost_chp1_pool,Cost_chp2_pool,Cost_p,Cost_r,Cost_seg,Name_pool,Name_r]=f_esm_dd(n,th_1,th_n,Change);
clear f_esm_dd
tic

NTC=repmat(Ntc,1,1,nh); % LATER: NTC should be imported as hourly network capacities from NordPool
n_dd=size(Hyd_infl_dd,1);   % Number of the days in the examined year
%-----------------------------------------------------------------------
% 2.2) External power market flows
% Power exchanges between NordPool and other exchanges
Exch_fi_ru=Market_exch(:,1);  %MWh/h
Exch_fi_ee=Market_exch(:,2);
Exch_se_de=Market_exch(:,3);
Exch_dk1_de=Market_exch(:,4);
Exch_dk2_de=Market_exch(:,5);
Exch_no_nl=Market_exch(:,6);
Exch_uk=Market_exch(:,7);

El_dem_exch=zeros(nh,n);         % Additional demand due to exchange with external markets (out of Nordpool) 
El_dem_exch(:,1)=Exch_fi_ru(th_1:th_n)+Exch_fi_ee(th_1:th_n);
El_dem_exch(:,2)=Exch_se_de(th_1:th_n);
El_dem_exch(:,3)=Exch_dk1_de(th_1:th_n);
El_dem_exch(:,4)=Exch_dk2_de(th_1:th_n);
El_dem_exch(:,5)=Exch_no_nl(th_1:th_n);
El_dem_exch(:,6)=Exch_uk(th_1:th_n);
%--------------------------------------------------------------------------------------
El_dem_pool=El_dem_res_tot+El_dem_exch;
El_dem_pool_uk=El_dem_res_tot(:,6)+El_dem_exch(:,6);
%-------------------------------------------------------------
% 2.3) Hydropower initialization 
np=length(Name_pool);
Hyd_res_d=zeros(n_dd,n);
Hyd_infl_d=zeros(n_dd,n);
Hyd_res_d(:,Region_wv)=Hyd_res_dd;
Hyd_infl_d(:,Region_wv)=Hyd_infl_dd;

% Rolling average of hydro inflow
for i=1:n_dd-7
    Hyd_infl_d(i,:)=mean(Hyd_infl_d(i:i+7,:),1);
end

El_pro_vre=El_pro_vre1+El_pro_vre2+El_pro_vre3+El_pro_vre4;

% 2.4) Building required matrixes  

bz_subindex = @(A,r,c) A(r,c);      % An anonymous function to index a matrix (value = bz_subindex(Matrix,row,column))

El_pro_hyd=zeros(nh,n);             % Actual hydropower production (MWh/h)
Hyd_reg_infl=zeros(n,nd);             % Regulated hydropower (total inflow for one day minus on-regulated hydro)
Stor_pre=zeros(n,nd);               % Storage residual from the previous day after bidding to the market
Stor_pre(:,1)=St_begin;             % Resrvoir level at the beginning of the first day 
Hyd_unaccepted=zeros(n,nd);                         % Residual hydro from production planning of the day before (if not sold in the market)
Hyd_storage=zeros(n,nd);                       % Unused hydro total

% El_dem_exch=zeros(n,nh);                       % The power demand offerd to the pool from each region (after domestic inflexible and non-market participant production)
Pow_act=zeros(n,np,nh);                            % The actual power production in each region-hour
Pow_mar=cell(n,nd);                            % The power production mix offerd to the pool from each region at each hour (not sorted, not guaranteed to be bought)
Cost_mar=cell(n,nd);                           % The cost matrix of power production mix for each hour-region presented to the pool (not sorted)
Pow_sort=cell(n,nh);                                % The sorted power production mix offerd to the pool from each region at each hour (sorted by their costs: supply curve)
Cost_sort=cell(n,nh);                               % The sorted cost matrix for each hour-region presented to the pool 

Exch_t=zeros(n,n,nh);                               % Creating matrixes for the results and post-processing
Sup_t=zeros(nh,n);
Area_p=zeros(nh,n);
Sys_p=zeros(nh,1);
Area_p_uk=zeros(nh,1);
Area_p_uk2=zeros(nh,1);

El_imp_crit=zeros(nh,n);                         %LATER should be defined how to supply deficit power after trade

% Caculation of Excess power
El_dem_res=El_dem_pool-El_pro_vre;      % Net residual demand: demand after fixed power, fixed import, and domestic VRE
El_exp_crit_dom=El_dem_res.*(El_dem_res<0);    %(-) sign         % DOMESTIC critical excess (without trade): demand after fixed power, fixed import, and domestic VRE

El_exp_crit_tech= El_exp_crit_dom+repmat((sum(Ntc,2))',nh,1);     %NOTICE: Ntc should be corrected if NTC is used                        % Critical excess assuming full-export possibilities: demand after fixed power, fixed import, and domestic VRE
El_exp_crit_tech(El_exp_crit_tech>0)=0;

El_exp_crit_mar=zeros(nh,n);     % Critical excess (after power market): demand after fixed power, fixed import, and domestic VRE, and trade with other countries

%_______________________________________________________________________________________________
% 3) Analysis 
% 3.1) Water value (Weekly)
N_week=zeros(52,1)+7;
f = @(k) repmat(Cost_hyd_corr(k,:), round(N_week(k)), 1);
Cost_hyd_corr_dd=cell2mat(arrayfun(f, (1:length(N_week))', 'UniformOutput', false));
Cost_hyd_corr_dd=[Cost_hyd_corr_dd;Cost_hyd_corr_dd(52*7-(n_dd-52*7)+1:end,:)];
Cost_hyd_corr_dd(isnan(Cost_hyd_corr_dd))=0;

El_pro_chp_fix=zeros(n,n);        % Share of CHP that CANNOT be reduced in excess VRE

dc=td_1-1;
disp('Ongoing simulations, please wait...!')
for d=1:nd              % d: count of day   
    h1=24*(d-1)+1;      % First hour of the analysis
    h2=24*d;            % Last hour of the analysis
    dc=dc+1;            % Counter of day numbers
    
    % 3.1.1) Hydro planning based on water value with a weekly horizon (using "hydroplan" function)
    
    Cost_hyd_wv=zeros(n,length(Ind_hyd));         % The number of cost segments in hydro water value curve (see function above)
    Hyd_pool_d=zeros(24,length(Ind_hyd),n);        % Hydro power (MWh/h) offered to the pool 
    
    for r=Region_wv      % Regions of hydropower water value
        Hyd_reg_infl(r,d)=Hyd_infl_d(dc,r)-sum(El_pro_vre4(h1:h2,r))/El_eff_hyd(r);   % IMPLEMENT EFFICIENCY IN OTHER INTERFACES
    
        [Hyd_mar_d, Cost_mod, stor_rem]=f_hydroplan(Cost_wv,Wv_seg,Cost_hyd_coff(:,r),Hyd_reg_infl(r,d),Stor_pre(r,d),St_min(r),Hyd_unaccepted(r,d),...
        El_dem_pool(h1:h2,r),El_pro_vre(h1:h2,r),...
        El_cap_hyd(r),St_hyd(r),El_eff_hyd(r));
        clear f_hydroplan
        Hyd_pool_d(:,:,r)=Hyd_mar_d;      % Hydro power sent to the mark for 24 h (with different water values)
        Cost_hyd_wv(r,:)=Cost_mod.*Cost_hyd_corr_dd(dc,r);
        if d<nd
           Stor_pre(r,d+1)=stor_rem;
        end
        Hyd_storage(r,d)=stor_rem;
    end
   %------------------------------------------------------------------------------------------------ 
    % 3.2) Construction of final power production portfolio to be sent to the day-ahead power market
    k=0;              % Counter of hours on every day
    I_sort=cell(n,24);
    
    for h=24*(d-1)+1:24*d                       % 24 h: hour of the day 
        k=k+1;
        % Name Convention: {'PV', 'Wind','Wind off', 'Riverhydro', 'Hydro', 'Nuclear',...
             %'CHP1 waste', 'CHP1 peat','CHP1 bio', 'CHP1 coal','CHP1 gas', 'CHP1 oil',...
             %'CHP2 waste', 'CHP2 peat','CHP2 bio', 'CHP2 coal','CHP2 gas', 'CHP2 oil',...
             %'Waste','Peat','Biomass','Coal','Gas','Oil','Reserve'};
        % 3.2.1) Power production matrix for bidding to the pool
      
        Pow_pool=[El_pro_vre1(h,:); El_pro_vre2(h,:); El_pro_vre3(h,:); El_pro_vre4(h,:) ;...
            squeeze(Hyd_pool_d(k,:,:)); El_pro_nuc(h,:);...
        El_chp1_mar(:,:,h);El_chp2_mar(:,:,h);El_ppchp_mar(:,:,h)];
      
        % 3.2.2) Cost matrix
        for r=1:n
         % PAY ATTENTION: the order of the cost elements in the vector below MUST conform to the order of power production mix in "Pow_pool" (no need for sorting based on the costs)          
            Cost_mar{r,h}=[Cost_p(1:4,r)' Cost_hyd_wv(r,:).*Cost_p(5,r)...
                Cost_p(6,r)' Cost_chp1_pool(:,r)' Cost_chp2_pool(:,r)' Cost_p(7:end,r)'];
         
            Pow_mar{r,h}=Pow_pool(:,r);
            Cost1=Cost_mar{r,h};
         
            II=find(fix(Pow_pool(:,r)));
            Power1=bz_subindex(Pow_pool(:,r),(fix(Pow_pool(:,r))>0),1);           % If some power production modes ar zero in that hour-region, they are excluded
            Cost1=Cost1(fix(Pow_pool(:,r))>0);
            [Cost_s, I_s]=sort(Cost1);           % Sorting the cost and power vectors ( for the creation of supply curve)
            Pow_sort{r,h}=Power1(I_s)';              % Constructing the power supply matrix for all the regions (cell)    
            Cost_sort{r,h}=Cost_s';
            I_sort{r,k}=II(I_s);                          % This is the matrix of indexs of NON-ZERO power production modes in Pow_pool sorted based on their cost
        end
        % The external power market modelled with Nord Pool through market coupling

        % Power price in the external market before trade with NordPool
        A1=cumsum(Pow_sort{6,h});
        [Sup_uk,Cost_uk ]=suppcurve_dis(A1, Cost_sort{6,h}, Cost_r(6),Cost_seg(6));
 

        i1=find(Sup_uk>=El_dem_pool(h,6),1);
        Area_p_uk(h)=Cost_uk(i1);
    end    
    %-----------------------------------------------------------------------------------
    % 3.3) Power market optimization
    % The matrixes of Demand, Cost and Power are sent ot the function
    % f_np_7R_dd so that the power market (including the network) will be
    % optimized. The outcome are area prices (Area_p), final production in
    % each region (Sup_tot), Optimal network flows (Exch_opt), the
    % information about the share of last producing unit, and system prices (sys_p).
    [Area_p0,Sup_tot, Exch_opt,Ind_last,Sup_last,Sys_p0]=f_np_dd(El_dem_pool(h1:h2,:)', Pow_sort(:,h1:h2),...
    Cost_sort(:,h1:h2),NTC(:,:,h1:h2),Cost_r,Cost_seg,Region_sys,peak_price);
    clear f_np_dd
    %----------------------------------------------------------------------------------- 
    % 3.4) Receieving the results and sorting them
    t=0;
    for h=24*(d-1)+1:24*d
        t=t+1;
        for r=1:n
              P_mode=zeros(size(Pow_mar{r,h},1),1);
              P=Pow_sort{r,h};
              P(Ind_last(r,t):end)=0; % Those power modes sent to the pool but not used become zero for actual production
              P(Ind_last(r,t))=Sup_last(r,t);
              P_mode(I_sort{r,t})=P;                         % Power production is set to the corresponding place in Pow_pool matrix
              Pow_act(r,:,h)=P_mode;
              El_pro_hyd(h,r)=sum(Pow_act(r,Ind_hyd,h))+El_pro_vre4(h,r);         
              El_exp_crit_mar(h,r)=min((Sup_tot(r,t)-El_pro_vre(h,r)),0);   % Critical excess power (with negative sign)
        end
    end
    %_______________________________________________________________________________________  
    % 4) General results
    Area_p(h1:h2,:)=Area_p0';
    Sys_p(h1:h2)=Sys_p0;
    Sup_t(h1:h2,:)=Sup_tot';          % Power that must be supplied internally after the market clearance and calculating power imports or EXPORTS
    Exch_t(:,:,h1:h2)=Exch_opt;                          
 
    for r=1:n
        if d<nd
            Hyd_unaccepted(r,d+1)=sum(sum(Hyd_pool_d(:,:,r),2))-sum(El_pro_hyd(h1:h2,r))/El_eff_hyd(r);
            Hyd_storage(r,d)=Hyd_storage(r,d)+Hyd_unaccepted(r,d+1);   % The content of storage from previous day + the part of hydro not accepted in the marker
        end
    end
   
    disp(['Day ' num2str(d) ' is complete!']) 
end
disp('Analysis for the whole period is successfully completed!')
toc
%_____________________________________________________________________________
% 5) Postprocessing and visulaization
Name_table=cell(length(Name_pool)+1,1);
Name_table(1:9)=Name_pool(1:9);Name_table{10}='(Total hydro)';Name_table(11:end)=Name_pool(10:end);
P_act_hyd=array2table(fix([Pow_act(:,1:9,h)'; sum(Pow_act(:,4:9,h),2)'; Pow_act(:,10:end,h)']),'RowNames',Name_table,'VariableNames',Name_r);

% 5.1) Loading historical data for comparison
A=load('Input_data.mat');
El_chp_his_fi=table2array(A.El_chp_fi);
El_hyd_his=table2array(A.El_hyd_pro);

Area_p_his=table2array(A.Area_p_hh);                            % Area prices (€/MWh) 
Sys_p_his=table2array(A.Sys_p_hh); 
%-------------------------------------------------
% 5.2) Power marker results
Pow_act=permute(Pow_act,[2 1 3]);
El_act_imp1=zeros(nh,n);
for i=1:n
El_act_imp1(:,i)=sum(Exch_t(:,i,:),1)-sum(Exch_t(i,:,:),2);
end
El_act_exp=-El_act_imp1.*(El_act_imp1<0);
El_act_imp=El_act_imp1.*(El_act_imp1>0);
El_act_pv=reshape(Pow_act(1,:,:),n,nh)';
El_act_wind=reshape(Pow_act(2,:,:),n,nh)';
El_act_windoff=reshape(Pow_act(3,:,:),n,nh)';
El_act_riverhyd=reshape(Pow_act(4,:,:),n,nh)';
El_act_hyddam=reshape(sum(Pow_act(5:5+length(Ind_hyd)-1,:,:)),n,nh)';
El_act_hydro=reshape(sum(Pow_act(4:5+length(Ind_hyd)-1,:,:)),n,nh)';
El_act_nuc=reshape(Pow_act(10,:,:),n,nh)';
El_act_chp1=reshape(sum(Pow_act(11:16,:,:)),n,nh)';
El_act_chp2=reshape(sum(Pow_act(17:22,:,:)),n,nh)';
El_act_chp=El_act_chp1+El_act_chp2;
% El_act_pp_waste=reshape(Pow_act(23,:,:),n,nh)';
% El_act_pp_peat=reshape(Pow_act(24,:,:),n,nh)';
% El_act_pp_bio=reshape(Pow_act(25,:,:),n,nh)';
% El_act_coal=reshape(Pow_act(26,:,:),n,nh)';
% El_act_gas=reshape(Pow_act(27,:,:),n,nh)';
% El_act_oil=reshape(Pow_act(28,:,:),n,nh)';
El_act_pp=reshape(sum(Pow_act(23:28,:,:)),n,nh)';
El_act_reserve=reshape(Pow_act(29,:,:),n,nh)';
El_act_ppchp1=min(Pow_act(23:28,:,:),El_ppchp2_mar);
El_act_ppchp=reshape(sum(El_act_ppchp1),n,nh)';
%------------------------------------------------------------------
% 5.3) Annual values (TWh/a)
% Power market
area_p=mean(Area_p);
sys_p=mean(Sys_p);

% District heating
dh_pro_chp=sum(Dh_pro_chp)./1e6; 
dh_pro_boil=sum(Dh_pro_boil)./1e6;
dh_pro_hob=sum(Dh_pro_hob)./1e6;

% Power production including CHP
el_act_pv=sum(El_act_pv)./1e6;
el_act_wind=sum(El_act_wind)./1e6;
el_act_windoff=sum(El_act_windoff)./1e6;
el_act_riverhyd=sum(El_act_riverhyd)./1e6;
el_act_hyd_dam=sum(El_act_hyddam)./1e6;
el_act_hyd=sum(El_act_hydro)./1e6;
el_act_nuc=sum(El_act_nuc)./1e6;

el_act_chp1=sum(El_act_chp1)./1e6;
el_act_chp2=sum(El_act_chp2)./1e6;
el_act_chp=el_act_chp1+el_act_chp2;

el_act_pp=sum(El_act_pp)./1e6;
el_act_reserve=sum(El_act_reserve)./1e6;

el_act_imp=sum(El_act_imp)./1e6;
el_act_exp=sum(El_act_exp)./1e6;

el_dem_tot=sum(El_dem_tot)./1e6;
el_imp_fix=sum(-El_dem_exch)./1e6;                  % Fixed or simulated import from countries outside Nordpool
el_pro_fix=sum(El_pro_fix_hh)./1e6;                 % Fixed domestic production

% Technology-based annual values for PP and CHP
El_act_chp1_a=sum(Pow_act(11:16,:,:),3)./1e6;
El_act_chp2_a=sum(Pow_act(17:22,:,:),3)./1e6;

El_act_chp_a=El_act_chp1_a+El_act_chp2_a;
El_act_pp_a=sum(Pow_act(23:28,:,:),3)./1e6;

%----------------------------------------------------------------------------------------------------
% 5.4) Fuel consumptions and emissions (TWh_fuel/a)
f_nuc=el_act_nuc./Eff_el(6,:);
F_pp_a=El_act_pp_a./Eff_el(7:12,:);

F_chp1_a=El_act_chp1_a./Eff_el_chp1;
F_chp2_a=El_act_chp2_a./Eff_el_chp2;

f_boil=dh_pro_boil./Eff_chp_boil;
f_hob=dh_pro_hob./Eff_hob;

F_boil_a=repmat(f_boil,size(F_boil_r,1),1).*F_boil_r;
F_hob_a=repmat(f_hob,size(F_hob_r,1),1).*F_hob_r;

%-------------------------------------------------------
% 5.5) Carbon emissions (million tonne/a)
Emm_fac=Emm_fac(7:12,:);
emm_pp=sum(F_pp_a.*Emm_fac,1); 
emm_chp1=sum(F_chp1_a.*Emm_fac,1);
emm_chp2=sum(F_chp2_a.*Emm_fac,1);
emm_chp=emm_chp1+emm_chp2;
emm_boil=sum(F_boil_a.*Emm_fac,1);
emm_chpdh=emm_chp+emm_boil;               % Emmision from CHP plants (including back-up boilers)

emm_hob=sum(F_hob_a.*Emm_fac,1);
emm_tot=emm_pp+emm_chp+emm_boil+emm_hob;
%________________________________________________________________________________________
% Dealing with excess power
el_exp_crit_dom=sum(El_exp_crit_dom)./1e6;     % TWh/a
el_exp_crit_tech=sum(El_exp_crit_tech)./1e6;
el_exp_crit_mar=sum(El_exp_crit_mar)./1e6;

% Reducing from CHP (Note: later it can be divided to each CHP category)
El_act_chp_vre=(El_act_chp+El_exp_crit_mar).*((El_act_chp+El_exp_crit_mar)>0);
El_exp_crit_chp=(El_act_chp+El_exp_crit_mar).*((El_act_chp+El_exp_crit_mar)<0);       % Critical excess after reducing power from CHP

el_act_chp_vre=sum(El_act_chp_vre)./1e6;        % Power from CHP after reducing critical VRE TWh/a             
dh_pro_chp_vre=el_act_chp_vre./Pth_chp;
dh_pro_boil_vre=dh_pro_boil+(dh_pro_chp-dh_pro_chp_vre);      % Boiler usage after critical excess

el_act_chp2_vre=(el_act_chp_vre-el_act_chp1).*((el_act_chp_vre-el_act_chp1)>0);
el_act_chp1_vre=el_act_chp_vre-el_act_chp2_vre;

F_chp1_vre=(repmat(el_act_chp1_vre,size(F_chp1_r,1),1).*F_chp1_r)./Eff_el_chp1;
F_chp2_vre=(repmat(el_act_chp2_vre,size(F_chp2_r,1),1).*F_chp2_r)./Eff_el_chp2;
f_boil_vre=dh_pro_boil_vre./Eff_chp_boil;
F_boil_vre=repmat(f_boil_vre,size(F_boil_r,1),1).*F_boil_r;

emm_chp1_vre=sum(F_chp1_vre.*Emm_fac,1);
emm_chp2_vre=sum(F_chp2_vre.*Emm_fac,1);
emm_chp_vre=emm_chp1_vre+emm_chp2_vre;
emm_boil_vre=sum(F_boil_vre.*Emm_fac,1);
emm_chpdh_vre=emm_chp_vre+emm_boil_vre;   % Emmision from CHP plants (including back-up boilers)
%-------------------------------------------------------------------------
% 5.4) saving the results in a mat file
save([save_name '.mat'],'Area_p','Sys_p','Exch_t','Sup_t','Pow_mar','Cost_mar','Pow_act','El_dem_pool','Pow_sort','Cost_sort')

%-------------------------------------------------------------------------
% 5.5) Visualization
figure
plot(Sys_p);hold on;
plot(Sys_p_his(th_1:th_n))
legend('Modeled prices','Historical data 2016')
ylim([0 1.2* max(Sys_p_his(th_1:th_n))])
set(gca, 'Xtick',0:24*7:24*td_n, 'XtickLabel', fix(td_1/7):td_n/7,'XGrid','on');

figure
plot(El_act_hydro(:,1));hold on;plot(El_hyd_his(th_1:th_n,1));hold on;

title( 'Hourly hydro production Finland (MWh/h)');legend('Simulation by Enerallt (FI)','Historical data (FI 2016)')
xlabel('Weeks in the examined period');ylim([0 1.2* max(El_hyd_his(th_1:th_n,1))])
set(gca, 'Xtick',0:24*7:24*td_n, 'XtickLabel', fix(td_1/7):td_n/7,'XGrid','on');

figure
plot(El_act_hydro(:,2));hold on;plot(El_hyd_his(th_1:th_n,2));hold on;

title( 'Hourly hydro production Sweden (MWh/h)');legend('Simulation by Enerallt SE','Historical data (SE 2016)')
xlabel('Weeks in the examined period');ylim([0 1.2* max(El_hyd_his(th_1:th_n,2))]);
set(gca, 'Xtick',0:24*7:24*td_n, 'XtickLabel', fix(td_1/7):td_n/7,'XGrid','on');

figure
plot(El_act_hydro(:,5));hold on;plot(El_hyd_his(th_1:th_n,3));hold on;

title( 'Hourly hydro production Norway (MWh/h)');legend('Simulation by Enerallt (NO)','Historical data (NO 2016)')
xlabel('Weeks in the examined period');ylim([0 1.2* max(El_hyd_his(th_1:th_n,3))])
set(gca, 'Xtick',0:24*7:24*td_n, 'XtickLabel', fix(td_1/7):td_n/7,'XGrid','on');


% figure
% plot(El_act_chp(:,1));hold on;plot(El_chp_his_fi(th_1:th_n));hold on;
% title( 'Hourly power production from CHP-DH Finland (MWh/h)');legend('Simulation by Enerallt','Historical data (2014)')
% xlabel('Hours in the examined period');ylim([0 1.2* max(El_chp_his_fi(th_1:th_n))])
% 
fig_r=1;
X=th_1:th_n;
figure

% subplot(2,2,1)
% plot(El_dem_tot(:,fig_r)+El_dem_exch(:,fig_r),'LineStyle','-');hold on;xlim([th_1 th_n]);
% area(X',[ El_act_nuc(:,fig_r) El_pro_fix_hh(th_1:th_n,fig_r) ],'LineStyle','none');
% xlabel('Examined hours');ylabel('MWh/h');title('Inflexible production');
% legend('Power demand','Nuclear','Industry');bz_reorderLegend([3,1,2]);
% 
% subplot(2,2,2)    
% plot(El_dem_tot(:,fig_r)+El_dem_exch(:,fig_r),'LineStyle','-');hold on;xlim([th_1 th_n]);
% area(X',[El_pro_vre4(:,fig_r) El_pro_hyd(:,fig_r)-El_pro_vre4(:,fig_r) El_pro_vre1(:,fig_r) El_pro_vre2(:,fig_r) El_pro_vre3(:,fig_r) ],'LineStyle','none');xlim([th_1 th_n])
% xlabel('Examined hours');ylabel('MWh/h');title('Power from RES plants (except bioenergy)');
% legend('Power demand','Riverhydro', 'Hydropower (storage)','Solar PV','Wind onshore','Wind-offshore');bz_reorderLegend([6,1,2,3,4,5]);
% xlim([th_1 th_n]);
% subplot(2,2,3)
% plot(El_dem_tot(:,fig_r)+El_dem_exch(:,fig_r),'LineStyle','-');hold on;xlim([th_1 th_n]);
% area(X',[El_chp_ind_hh(th_1:th_n,fig_r) El_act_chp1(:,fig_r) El_act_chp2(:,fig_r) El_act_pp(:,fig_r)],'LineStyle','none');
% xlabel('Examined hours');ylabel('MWh/h');title('Thermal power plants');
% legend('Power demand','Industrial CHP','Local CHP','Central CHP','Power plant');bz_reorderLegend([5,1,2,3,4]);
% 
% subplot(2,2,4)
plot(X',El_dem_tot(:,fig_r)+El_dem_exch(:,fig_r),'LineWidth',0.3,'Color','black');hold on;xlim([th_1 th_n]);           % The name should be adjusted based on each desirable production mode
area(X',[El_act_nuc(:,fig_r) El_pro_hyd(:,fig_r),...
      El_pro_fix_hh(th_1:th_n,fig_r) El_act_chp1(:,fig_r)+El_act_chp2(:,fig_r),...
     El_pro_vre2(:,fig_r)+El_pro_vre3(:,fig_r) El_pro_vre1(:,fig_r),...
     El_act_pp(:,fig_r) El_act_imp(:,fig_r)],'LineStyle','none');
xlabel('Examined hours');ylabel('Power supply/demand (MWh/h)');title(['Power supply mix in region No. ' num2str(fig_r) ]);
legend('Power demand','Nuclear','Hydropower (total)','Industry (CHP)','CHP-DH','Wind','Solar PV', 'Condensing plants', 'Power imports');
chh = get(gca,'Children');
bz_reorderLegend([9,1,2,3,4,5,6,7,8],chh);hold on;
save(save_name,'-v7.3');
