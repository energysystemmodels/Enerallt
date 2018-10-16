% Input data file
%-------------------------------------------------------------------------------------------------------
% P: Power capacity
% T: Thermal capacity
% D: District heating 
% E: Electricity 
% d: demand
% p: production
%--------------------------------------------------------------------------------------------------
function [El_dem_hh,P_pow,Ramp,Eff_el,Emm_fac,El_chp_ind_hh,El_pro_fix_hh,El_pro_nuc_r,...
 Dd_hob_hh,Dd_chp_ind_hh,Dd_chp1_hh,Dd_chp2_hh,Chp_ind_hh,P_chp,P_chp1,Eff_el_chp1,Eff_el_chp2,Eff_th_chp1,Eff_th_chp2,Eff_chp_boil,Eff_hob,...
 Pth_chp,F_hob_r,F_chp1_r,F_chp2_r,F_boil_r,P_hp_chp1,P_hp_chp2,Cop_hp_chp1,Cop_hp_chp2,T_chp_boil,Chp_boil_fix,...
 St_chp1,St_chp2,Dp_sol_hob,Dp_sol_chp1,Dp_sol_chp2,Dp_sol_hob_r,Dp_sol_chp1_r,Dp_sol_chp2_r,St_sol_hob,St_sol_chp1,St_sol_chp2,Loss_sol_hob,Loss_sol_chp1,Loss_sol_chp2,... 
 Cf_pv,Cf_wind,Cf_windoff,Pv_hh,Wind_hh,Hyd_res_hh,Hyd_infl_hh,El_rivhyd_pro,Nuc_hh,Hyd_res_dd,Hyd_infl_dd,St_hyd,Cost_hyd_corr,Wv_seg,Cost_wv,Cost_hyd_coff,St_begin,St_min,...
 Name_r,Name_p,Name_chp,Ntc,Market_exch,Cost_p,Cost_chp1_mc,Cost_chp2_mc,Cost_ppchp2_mc,Cost_r,Cost_seg]=f_inputdata_r(n_tot)

Name_r= {'FI' 'SE' 'DK1' 'DK2' 'NO' 'UK'};  % Name of regions (used for table row names)
A=load('Input_data.mat');
%__________________________________________________________________________________________________________________
% 1) Network information and energy demand in each region

% 1.1) Net transmisson capacity between areas (bidirectional)
% Columns: to FI	to SE	to DK1 to DK2	to NO	to DE	to EE	to RU	to PO	to NE	to UK
% Rows: from FI, from SE, from DK1, from DK2, from NO, from DE, from EE, from RU, from PO, from NE, from UK

Ntc=[0    2300  0     0     90    0   ;    
    2700  0     680   1300  3995  0   ;
    0     740   0     590   1632  0   ;
    0     1700  600   0     0     0   ;
    90    3695  1632  0     0     1400;
    0     0     0     0     1400  0   ];
n_region=6;
%------------------------------------------------------------------------------
% 1.2) Electricity and district heating (DH) demand 
% 1.2.1) Electricity demand

El_dem_hh=table2array(A.El_dem);                                          % Hourly electricity demand (MWh/h) 2014
El_dem_uk=load('uk_eldem_2014.txt');

El_dem_hh=[El_dem_hh El_dem_uk];

nhp=size(El_dem_hh,1);                                                    % Number of hours in the examined period (year)
El_chp_ind_hh=table2array(A.El_pro_chp_ind);                              % Hourly electricity production from CHP plants in Finland (MWh/h) (will be used for industry fix production pattern)

Chp_ind_hh=El_chp_ind_hh./max(El_chp_ind_hh);                             % Relative hourly distribution of industrial CHP 

Ep_chp_ind=[9.28  5.8    (6028/(6028+990))*1.9     (990/(6028+990))*1.9    0  0];     % Annual electric production fom industry or autoproducers in CHP mode (TWh/a)

El_chp_ind_hh=repmat(Ep_chp_ind.*1e6,nhp,1).*repmat(Chp_ind_hh./sum(Chp_ind_hh),1,n_region);              % Hourly electric production fom industry or autoproducers in CHP mode (MWh/h)                        

Ep_fixed_other=[0   0    0     0     0  0];                          % Annual fixed electricity production (OTHER THAN INDUSTRY CHP)(TWh/a)
Ep_fixed_other_hh=repmat(Ep_fixed_other.*1e6/nhp,nhp,1);           % Hourly fixed production (here is constant throughout the year) (MWh/h)

% Ep_fixed_tot=Ep_chp_ind+Ep_fixed_other;                          % Total annual electricity production from NON-ADJUSTIBLE and autoproducer plants (TWh/a) 
El_pro_fix_hh=El_chp_ind_hh+Ep_fixed_other_hh;
%-----------------------------------------------------------------------------------------------------
% 1.2.1.1) Electricity prices (for validation and comparison)
Area_p_his_hh=table2array(A.Area_p_hh);                            % Area prices (€/MWh) 2014
Sys_p_his_hh=table2array(A.Sys_p_hh);                              % System price (€/MWh) 2014

%------------------------------------------------------------------------------------------------------------------------
% 1.2.2) District heating (DH) 
Dd=[34.7    51.7    (6028/(6028+990))*37.53     (990/(6028+990))*37.53     4.6   11.7];                     % Annual district heating production (TWh_th/a), Danish areas based on CHP capacity share
Dd_loss=[34.7-31.6  '?'  (6028/(6028+990))*7.4     (990/(6028+990))*7.4      '?'];                      % Annual DH distribution losses (TWh/a)
Dd_chp_all=[0.724*Dd(1)   0.4*57.2   (6028/(6028+990))*27.3     (990/(6028+990))*27.3  0.5*4.6  11.7];             % Annual DH demand that can be met by CHP (otherwise heat-only boilers) (TWh_th/a)
Dd_ind1=[0.1*Dd(1)  0.1*Dd(2)   (6028/(6028+990))*(4.9+1.5)     (990/(6028+990))*(4.9+1.5)      0.1*4.6  0];    % Annual DH demand of industry or autoproducers for DK (Twh/a)
Dd_chp_ind1=[0.1*Dd_chp_all(1)  0.1*Dd_chp_all(2)   (6028/(6028+990))*(4.9)     (990/(6028+990))*(4.9)      0.1*Dd_chp_all(5) 0];  % Annual DH demand of CHP industry or CHP autoproducers for DK (Twh/a)

Dd_chp_fixed_other=[ 0 0 0 0 0 0];                                   % Other fixed DH production in connection with CHP (e.g., from waste plants)

Dd_chp_ind=Dd_chp_ind1+Dd_chp_fixed_other;                         % Annual non-adjustible CHP

Dd_hob=Dd-Dd_chp_all;                                              % Annual DH demand that can be met by heat-only boilers (TWh_th/a)
Dd_chpdh=Dd_chp_all-Dd_chp_ind;                                    % Annual DH demand that can be met by CHP (excluding CHP-autonomous) (TWh_th/a)

Dd_chp1_r=[0.41   0.72   0.24   1     0.59  1];                      % Maximum share of DH demand that can be met by local CHP (based on capacity of small CHP relative to total P_chp), from CHP-DH, not including industry)

% Hourly DH demand modelling
Temp_hh=table2array(A.Temp);
Temp_uk=load('uk_temp_2014.txt');

Temp_hh=[Temp_hh Temp_uk];

El_chp_fi=table2array(A.El_chp_fi);

T_min=[15.5    15.5     15.5    15.5    15.5  15.5];                  % Min outdoor temperature (C) for space heating start running 
Hw_r= [0.3     0.3      0.3     0.3     0.3   0.3];                  % Share of hot water demand from annual heat demand

Hdh_hh=(repmat(T_min,size(Temp_hh,1),1)-Temp_hh).*((repmat(T_min,size(Temp_hh,1),1)-Temp_hh)>0)+...
    repmat(sum((repmat(T_min,size(Temp_hh,1),1)-Temp_hh).*((repmat(T_min,size(Temp_hh,1),1)-Temp_hh)>0),1).*Hw_r./size(Temp_hh,1),size(Temp_hh,1),1);       % Hourly degree hours
Hdh_r_hh=Hdh_hh./repmat(max(Hdh_hh),size(Hdh_hh,1),1);         % Relative HDH

Dd_hh=Hdh_r_hh.*(1e6*repmat(Dd,size(Hdh_hh,1),1)./repmat(sum(Hdh_r_hh,1),size(Hdh_hh,1),1));    % Hourly DH demand (MWh_th/h)

Dd_chp_ind_hh=repmat(Dd_chp_ind.*1e6,nhp,1).*repmat(Chp_ind_hh./sum(Chp_ind_hh),1,n_region);     % Hourly DH production of CHP-industry or CHP autoproducers (Mwh_th/h), (NON_ADJUSTIBLE CHP)

Dd_chp_hh=(Dd_hh-Dd_chp_ind_hh).*repmat(Dd_chpdh./(Dd_hob+Dd_chpdh),nhp,1);    % Hourly DH demand for (non-industrial) CHP plants (MWh_th/h)
Dd_hob_hh=fix(Dd_hh-Dd_chp_ind_hh-Dd_chp_hh);                                  % Hourly DH demand for heat-only plants (MWh_th/h)

%%%%%%%%% FINLAND
Dd_chp_hh(:,1)=El_chp_fi./(max(El_chp_fi)).*((Dd_chpdh(1)*1e6))/sum(El_chp_fi./(max(El_chp_fi)));

% Solar DH
Dp_sol_hob= [0   0   0   0    0  0];                            % Annual solar thermal DH production, connected to heat-only DH (TWh_th/a)
Dp_sol_chp1=[0   0   0   0    0  0];                            % Annual solar thermal DH production, connected to small CHP-DH (TWh_th/a)
Dp_sol_chp2=[0   0   0   0    0  0];                            % Annual solar thermal DH production, connected to large CHP-DH (TWh_th/a)


Dp_sol_hob_r=[0.8   0.8   0.8   0.8    0.8  0.8];                  % Maximum share of DH connected to heat-only boilers that can be met by solar thermal DH production (%)
Dp_sol_chp1_r=[0.8   0.8   0.8   0.8    0.8  0.8];                  % Maximum share of DH connected to small CHP-DH that can be met by solar thermal DH production (%)
Dp_sol_chp2_r=[0.8   0.8   0.8   0.8    0.8  0.8];                  % Maximum share of DH connected to large CHP-DH that can be met by solar thermal DH production (%)

Loss_sol_hob=[0.1   0.1   0.1   0.1    0.1  0.1];                   % Loss solar production connected to heat-only boilers (storage and etc) (%)
Loss_sol_chp1=[0.1   0.1   0.1   0.1    0.1  0.1];
Loss_sol_chp2=[0.1   0.1   0.1   0.1    0.1  0.1];
%________________________________________________________________________________________________________________________
% 2) Technical characteristics and availability of production plants
% 2.1) Installed production capacity
% 2.1.1) Only-power production capacity 

Name_p={'PV', 'wind on', 'wind off', 'riverhydro', 'hydro', 'nuclear',...
'cond. waste', 'cond. peat', 'cond. biomass', 'cond. coal', 'cond. gas', 'cond. oil', 'Reserve'};   
   
P_pow=[   10          45             421         180        0              5400            % 1. solar PV
          496         3550           2966        608        820            12440-4494      % 2. Onshore wind
          0           1925           843         428        0              4494      % 3. Offshore wind
          924         5825           0           0          6255           1570      % 4. River hydro (unregulated) % SE 5825 MW and NO 6255 MW
          2640-300        16650-5825     4           5          30136-6255 0      % 5. Hydro   (estimations for riverhydro is deducted)
          2752        8800           0           0          0              9940      % 6. Nuclear
          0           7              27          0          3              0      % 7. Waste
          120         5              0           0          0              0      % 8. Peat
          0           231            19          0          0              2120      % 9. Biomass thermal
          1047        64             1306        34         8              19900     % 10. Coal thermal
          128         238+410        69          11         877            31990+1790      % 11. Gas thermal (gas +cc)
          1560        2364           808         201        96             1370+1180      % 12. Oil thermal
          1600        1600           1600        1600       1600         1600 ];     % 13. Reserve plants (for emergency cases)
 
%-------------------------------------------------------------------------------------------------------------------------
% 2.1.2) CHP installed capacity in each area (MW_e), 

Name_chp={'CHP waste', 'CHP peat', 'CHP biomass', 'CHP coal', 'CHP gas', 'CHP oil'};

P_chp=[ 20+7            286        230         30          84        0.10*2066        % Waste      % For finland: CHP DH + CHP industry
        997+550         141        0           0           0           0      % Peat
        155+1610        2305       285         88          0           0.11*2066      % Biomass
        1596+20       395        3260        0           4           0.04*2066      % Coal
        1220+360        208        2251        821         547         0.67*2066      % Gas
        21+57          874        2           52          12        0.09*2066];      % Oil
%--------------------------------------------------------------------------------------------------------
% 2.1.3) Ramping constraints of power and CHP plants (as percentage of
% production in previous time step (hour)) for VRE, a large number of 100
Ramp=[1e2 1e2 1e2 1e2 1e2 1e2             % 'PV'
      1e2 1e2 1e2 1e2 1e2 1e2             %'Wind'
      1e2 1e2 1e2 1e2 1e2 1e2             %'Wind off'
      1e2 1e2 1e2 1e2 1e2 1e2             %'Riverhydro' 
      0.9 0.9 0.9 0.9 0.9 0.9             % 'Hydro'
      0.5 0.5 0.5 0.5 0.5 0.5             % 'Nuclear'
      0.5 0.5 0.5 0.5 0.5 0.5             %'CHP1 waste'
      0.6 0.6 0.6 0.6 0.6 0.6             %'CHP1 peat'
      0.9 0.9 0.9 0.9 0.9 0.9             % 'CHP1 bio'
      0.5 0.5 0.5 0.5 0.5 0.5             %'CHP1 coal'
      0.7 0.7 0.7 0.7 0.7 0.7             %'CHP1 gas'
      0.9 0.9 0.9 0.9 0.9 0.9              %'CHP1 oil'
      0.5 0.5 0.5 0.5 0.5 0.5             %'CHP2 waste'
      0.6 0.6 0.6 0.6 0.6 0.6             %'CHP2 peat'
      0.9 0.9 0.9 0.9 0.9 0.9             % 'CHP2 bio'
      0.5 0.5 0.5 0.5 0.5 0.5             %'CHP2 coal'
      0.7 0.7 0.7 0.7 0.7 0.7             %'CHP2 gas'
      0.9 0.9 0.9 0.9 0.9 0.9              %'CHP2 oil'
      0.9 0.9 0.9 0.9 0.9 0.9              %'Waste'
      0.6 0.6 0.6 0.6 0.6 0.6              %'Peat'
      0.5 0.5 0.5 0.5 0.5 0.5             %'Biomass'
      0.5 0.5 0.5 0.5 0.5 0.5               %'Coal'
      0.9 0.9 0.9 0.9 0.9 0.9              %'Gas'
      0.9 0.9 0.9 0.9 0.9 0.9              %'Oil'
      0.9 0.9 0.9 0.9 0.9 0.9];             %'Reserve'};
Ramp = Ramp.*1;    
    
% 2.1.3) Back-up boilers to supplement CHP plants when the DH demand exceeds CHP capacitiy [capcatires are arbitrary here)

T_chp_boil=[13000   20000   5000   6000   5000 5000];         % Thermal capacity of the back-up boilers (MW_th)
Chp_boil_fix=[0   0   0   0   0     0];                           % Fixed boiler ratio in connection with CHP (%)

F_chp_boil=[0     0     0     0     0    0            % Biomass 
            0     0     0     0     0    0            %  Coal 
            0.5   0.5   0.5   0.5   0.5   0.5            %  Gas 
            0.5   0.5   0.5   0.5   0.5  0.5 ];         % Oil

% 2.1.4) Capacity of local CHP (lower than 100 MW from Platts)
% CHP1: small and local CHP (just based on heat demand)
% CHP2: Large and central CHP (equipped with condenser)

P_chp1=[ 20+7          286        230         30          84      0.10*2066-1           % Waste      % For finland: CHP DH + CHP industry
        444+364         21         0           0           0         0        % Peat
        238+1129        1952       285         88          0          0.11*2066-1        % Biomass
        0+20         181        0           0           4         0.04*2066-1         % Coal
        195+232         208        934         812         287       0.67*2066-1         % Gas
        21+57          404        2           52          12       0.09*2066-1  ];      % Oil
%----------------------------------------------------------------------------------------
% 2.1.5) Large heat pumps and storage capacity (connected to DH)

P_hp_chp1=[0      0      0      0      0      0];            % Electric capacity of LHP in connection with small CHP-DH (MWe) 
P_hp_chp2=[0      0      0      0      0      0]; 

St_chp1=[0      0      0      0      0     0 ];               % Thermal storage capacity of small CHP-DH (MWh_th)
St_chp2=[0      0      0      0      0     0 ];               % Thermal storage capacity of large CHP-DH (MWh_th)                     

%----------------------------------------------------------------------------------------
% 2.2) Storage technologies and capacities 
% 2.2.1) Hydro storage
Hyd_reserve=table2array(A.Hyd_reserve);                      % Weekly hydro reservoir level in the Nordics (GWh) (NOTICE the unit)
Hyd_inflow=table2array(A.Hyd_inflow);                      % Weekly hydro inflow (GWh) (NOTICE the unit)

Hyd_res=Hyd_reserve.*1000;                    %MWh
%Hyd_res=[Hyd_res(:,1:2) zeros(52,2) Hyd_res(:,3) zeros(52,1)];   % NOTICE: added to make the number of columns equal to number of region
Hyd_infl=Hyd_inflow.*1000;                    %MWh
%Hyd_infl=[Hyd_infl(:,1:2) zeros(52,2) Hyd_res(:,3) zeros(52,1)];
% Constructing hourly hydro reserve from weekly (equal distribution)
N_week=zeros(52,1)+7*24;
f = @(k) repmat(Hyd_res(k,:)./(7*24), round(N_week(k)), 1);
Hyd_res_hh = cell2mat(arrayfun(f, (1:length(N_week))', 'UniformOutput', false));
Hyd_res_hh=[Hyd_res_hh;Hyd_res_hh(8736-23:8736,:)];
Hyd_res_hh(isnan(Hyd_res_hh))=0;

f = @(k) repmat(Hyd_infl(k,:)./(7*24), round(N_week(k)), 1);
Hyd_infl_hh = cell2mat(arrayfun(f, (1:length(N_week))', 'UniformOutput', false));
Hyd_infl_hh=[Hyd_infl_hh;Hyd_infl_hh(8736-23:8736,:)];
Hyd_infl_hh(isnan(Hyd_infl_hh))=0;

N_week=zeros(52,1)+7;
f = @(k) repmat(Hyd_res(k,:)./7, round(N_week(k)), 1);
Hyd_res_dd=cell2mat(arrayfun(f, (1:length(N_week))', 'UniformOutput', false));
Hyd_res_dd=[Hyd_res_dd;Hyd_res_dd(364,:)];
Hyd_res_dd(isnan(Hyd_res_dd))=0;

f = @(k) repmat(Hyd_infl(k,:)./7, round(N_week(k)), 1);
Hyd_infl_dd=cell2mat(arrayfun(f, (1:length(N_week))', 'UniformOutput', false));
Hyd_infl_dd=[Hyd_infl_dd;Hyd_infl_dd(364,:)];
Hyd_infl_dd(isnan(Hyd_infl_dd))=0;

St_hyd=1e6.*[4   33.7   0   0    85 4];             % Hydro storage MWh
St_begin=1e3*[3617 22406  0   0   55198 0];         % The reservoir level (MWh) at the beginning of the analysis
St_min=1e3*[1915   8030   0   0  26876 0];          % The minimum storage level (MWh) throughout year (from historical data)
% River Hydro storage
El_rivhyd_pro=table2array(A.El_rivhyd_pro);                      % Weekly hydro reservoir level in the Nordics (GWh) (NOTICE the unit)

%-----------------------------------------------------------------------------------------
% 2.2.2) Storage for olar thermal DH
St_sol_hob= [0   0   0   0    0 0];                 % Thermal storage capacity of solar thermal DH production, connected to heat-only DH (MWh_th)
St_sol_chp1=[0   0   0   0    0 0];                 % Thermal storage capacity of solar thermal DH production, connected to small CHP-DH (MWh_th)
St_sol_chp2=[0   0   0   0    0 0];                 % Thermal storage capacity of solar thermal DH production, connected to large CHP-DH (MWh_th)

% 2.2.3) Storage for CHP-DH
St_chp1=[7000   3000   10000*(6028/(6028+990))   10000*(990/(6028+990))    0 0];                            % Thermal storage capacity small CHP-DH (MWh_th)
St_chp2=[10000   4000   50000*(6028/(6028+990))   50000*(990/(6028+990))   0 0];                            % Thermal storage capacity large CHP-DH (MWh_th)


%-----------------------------------------------------------------------------------------------------------------------
% 2.3) Efficiencies and COPs
% 2.3.1) Electric efficiency of the power-only plants  (for VRE, the inefficiencies can be reflected either in availability factors or here)

Eff_el=        [1    1    1    1    1     1              % 1. solar PV
                1    1    1    1    1     1              % 2. Onshore wind
                1    1    1    1    1     1              % 3. Offshore wind
                0.88   0.88   0.85   0.85   0.88  0.88                  % 4. River hydro (unregulated)
                0.88   0.88   0.85   0.85   0.88   0.88               % 5. Hydro
                0.33   0.33   0.33   0.33   0.33    0.33              % 6. Nuclear
                0.23   0.23   0.23   0.23   0.23    0.23              % 7. Waste
                0.28   0.28   0.28   0.28   0.28    0.28              % 8. Peat
                0.35   0.35   0.35   0.35   0.35    0.35              % 9. Biomass thermal
                0.41   0.41   0.41   0.41   0.41     0.41             % 10. Coal thermal
                0.45   0.45   0.45   0.45   0.45     0.45             % 11. Gas thermal
                0.41   0.41   0.41   0.41   0.41     0.41             % 12. Oil thermal
                0.41   0.41   0.41   0.41   0.41    0.41];                % 13. Reserve plants (for emergency cases)
                
% 2.3.2) Electric efficiency of CHP plants  
Eff_el_chp1=    [0    0    0    0    0     0              % 1. solar PV
                 0    0    0    0    0     0              % 2. Onshore wind
                 0    0    0    0    0     0              % 3. Offshore wind
                 0    0    0    0    0     0              % 4. River hydro (unregulated)
                 0    0    0    0    0     0              % 5. Hydro
                 0    0    0    0    0     0             % 6. Nuclear
                 0.18   0.18   0.18  0.18    0.18  0.18        % 7. Waste
                 0.2    0.2    0.2   0.2     0.2    0.2               % 8. Peat
                 0.24   0.24   0.24   0.24   0.24   0.24               % 9. Biomass thermal
                 0.25   0.25   0.25   0.25   0.25    0.25              % 10. Coal thermal
                 0.28   0.28   0.28   0.28   0.28    0.28              % 11. Gas thermal
                 0.28   0.28   0.28   0.28   0.28    0.28              % 12. Oil thermal
                 0      0      0      0      0   0];               % 13. Reserve plants (for emergency cases)

Eff_el_chp2=    [0    0    0    0    0      0             % 1. solar PV
                 0    0    0    0    0      0             % 2. Onshore wind
                 0    0    0    0    0      0             % 3. Offshore wind
                 0    0    0    0    0      0             % 4. River hydro (unregulated)
                 0    0    0    0    0      0             % 5. Hydro
                 0    0    0    0    0       0            % 6. Nuclear
                 0.2    0.2    0.2    0.2    0.2   0.2      % 7. Waste
                 0.22   0.22   0.22   0.22   0.22   0.22               % 8. Peat
                 0.25   0.25   0.25   0.25   0.25   0.25               % 9. Biomass thermal
                 0.28   0.28   0.28   0.28   0.28    0.28              % 10. Coal thermal
                 0.30   0.30   0.30   0.30   0.30    0.3              % 11. Gas thermal
                 0.30   0.30   0.30   0.30   0.30    0.3              % 12. Oil thermal
                 0      0      0      0      0   0 ];               % 13. Reserve plants (for emergency cases)

Pth_chp=[0.65 0.5 0.5 0.5 0.5 0.5];                                      % Power ot heat in national level (based on Finnish DH data Energiateollisuus)
    
%----------------------------------------------------------------------------------------------------
% 2.3.3) COP of large heat pumps (connected to DH)

Cop_hp_chp1=[3      3      3      3      3 3];            % COP of LHP in connection with small CHP-DH (MWe) 
Cop_hp_chp2=[3      3      3      3      3 3]; 

%-----------------------------------------------------------------------------------
% 2.2.4) Thermal efficiency of CHP plants  
Eff_th_chp1=    [0    0    0    0    0     0             % 1. solar PV
                 0    0    0    0    0     0             % 2. Onshore wind
                 0    0    0    0    0     0              % 3. Offshore wind
                 0    0    0    0    0     0             % 4. River hydro (unregulated)
                 0    0    0    0    0     0              % 5. Hydro
                 0    0    0    0    0     0              % 6. Nuclear
                 0.4   0.4   0.4   0.4   0.4    0.4              % 7. Waste
                 0.4   0.4   0.4   0.4   0.4    0.4              % 8. Peat
                 0.5   0.5   0.5   0.5   0.5    0.5               % 9. Biomass thermal
                 0.5   0.5   0.5   0.5   0.5    0.5               % 10. Coal thermal
                 0.5   0.5   0.5   0.5   0.5    0.5               % 11. Gas thermal
                 0.5   0.5   0.5   0.5   0.5    0.5               % 12. Oil thermal
                 0     0     0     0     0 0];                     % 13. Reserve plants (for emergency cases

Eff_th_chp2=    [0    0    0    0    0     0             % 1. solar PV
                 0    0    0    0    0     0             % 2. Onshore wind
                 0    0    0    0    0     0              % 3. Offshore wind
                 0    0    0    0    0     0             % 4. River hydro (unregulated)
                 0    0    0    0    0     0              % 5. Hydro
                 0    0    0    0    0     0              % 6. Nuclear
                 0.4   0.4   0.4   0.4   0.4    0.4              % 7. Waste
                 0.4   0.4   0.4   0.4   0.4    0.4              % 8. Peat
                 0.5   0.5   0.5   0.5   0.5    0.5               % 9. Biomass thermal
                 0.5   0.5   0.5   0.5   0.5    0.5               % 10. Coal thermal
                 0.5   0.5   0.5   0.5   0.5    0.5               % 11. Gas thermal
                 0.5   0.5   0.5   0.5   0.5    0.5               % 12. Oil thermal
                 0     0     0     0     0 0];                     % 13. Reserve plants (for emergency cases

Eff_chp_boil=[0.88   0.88   0.88   0.88   0.88 0.88];         % Thermal efficiency of the back-up boilers (%)
Eff_hob=     [0.82   0.82   0.85   0.85   0.88 0.88];         % Thermal efficiency of the DH heat-only boilers (%)
%---------------------------------------------------------------------------------------------------------------------
% 2.4) Availability of variable renewable energy and other power production plants

Cf_pv=    [850    900    1000    1000    850  950];    % PV production kWh/kW per year
Cf_wind=   [0.24   0.23   0.27   0.28   0.23  0.30];    % Onshore wind capacity factors for simulation (slightly different from annual full-load factor!)
Cf_windoff=[0.27   0.26   0.28   0.30   0.28   0.34];    % Offshore wind capacity factors 

% Cf_pv=Cf_pv1/8760;
% Solar PV
Pv_hh=table2array(A.Solar);          % Hourly PV irradiation or production 2014
Pv_uk=load('uk_pv_2014.txt');
Pv_hh=[Pv_hh Pv_uk];
Pv_hh=Pv_hh./repmat(max(Pv_hh),nhp,1);                              % Relative

% Wind
Wind_hh=table2array(A.Wind);        % Hourly wind production or speed 2014
Wind_uk=load('uk_wind_2014.txt');
Wind_hh=[Wind_hh Wind_uk];

Wind_hh=Wind_hh./repmat(max(Wind_hh),nhp,1);                          % Relative
%-----------------------------------------------------------------------------
% Nuclear availability
Nuc_hh=table2array(A.Nuc_ava);          % Hourly nuclear availability 2014
Nuc_uk=load('uk_nuc_2014.txt');


Nuc_hh=[Nuc_hh repmat(Nuc_hh(:,1),1,n_region-3) Nuc_uk];                  % Just in case nuclear for other regions!
Nuc_hh=Nuc_hh./repmat(max(Nuc_hh),nhp,1);                          % Relative

El_pro_nuc_r=[0.3  0.3  0.3  0.3   0.3 0.3];                          % Share of nuclear power that cannot be lowered at all
%------------------------------------------------------------------------------------------------------------------------ 
% 2.5) Carbon emission factors (kg/kWh_fuel)  
Emm_fac=       [0    0    0    0    0   0                % 1. solar PV
                0    0    0    0    0    0               % 2. Onshore wind
                0    0    0    0    0   0                % 3. Offshore wind
                0    0    0    0    0   0                 % 4. River hydro (unregulated)
                0    0    0    0    0   0                % 5. Hydro
                0    0    0    0    0   0                 % 6. Nuclear
                0    0    0    0    0    0                % 7. Waste
                0.39    0.39    0.39    0.39    0.39    0.39                % 8. Peat
                0       0       0       0       0        0               % 9. Biomass thermal
                0.34    0.34    0.34    0.34    0.34     0.34               % 10. Coal thermal
                0.21    0.21    0.21    0.21    0.21      0.21              % 11. Gas thermal
                0.26    0.26    0.26    0.26    0.26     0.26               % 12. Oil thermal
                0.21    0.21    0.21    0.21    0.21   0.21];                  % 13. Reserve plants (for emergency cases)

%__________________________________________________________________________________________________________________________
% 3) Fuel used in energy production

F_hob=[0.12	  0.09	0.08	0.08	0.71    0   % waste
       0.11	  0.04	0.00	0.00	0.00    0   % peat
       0.37	  0.82	0.50	0.50	0.24    0.1   % Biomass
       0.03	  0.00	0.00	0.00	0.00     0.1  % Coal
       0.25	  0.02	0.38	0.38	0.04     0.7  % Gas
       0.09	  0.04	0.04	0.04	0.01     0.1];     % Oil
   
F_hob_r=F_hob./repmat(sum(F_hob,1),size(F_hob,1),1);             % Ratios

F_chp1=[0.03	0.38	0.21	0.21	0.71       % waste       % Small CHP
        0.14	0.03	0.00	0.00	0.00       % peat
        0.31	0.46	0.16	0.16	0.24       % Biomass
        0.25	0.07	0.01	0.01	0.00       % Coal
        0.22	0.05	0.54	0.54	0.04       % Gas
        0.02	0.01	0.08	0.08	0.01];     % Oil

F_chp1_r=P_chp1./repmat(sum(P_chp1,1),size(P_chp1,1),1);
    
F_chp2=[0.03	0.38	0.00	0.00	0.71       % waste       % Large CHP
        0.14	0.03	0.00	0.00	0.00       % peat
        0.31	0.46	0.18	0.18	0.24       % Biomass
        0.25	0.07	0.72	0.72	0.00       % Coal
        0.22	0.05	0.09	0.09	0.04       % Gas
        0.02	0.01	0.01	0.01	0.01];     % Oil
    
%F_chp2_r=F_chp2./repmat(sum(F_chp2,1),size(F_chp2,1),1);

F_boil=[0.00	0.00	0.00	0.00	0.00   0    % waste       % Large CHP
        0.00	0.00	0.00	0.00	0.00   0    % peat
        0.00	0.00	0.00	0.00	0.00    0   % Biomass
        0.00	0.00	0.00	0.00	0.00    0   % Coal
        0.50	0.50	0.50	0.50	0.50   0.5    % Gas
        0.50	0.50	0.50	0.50	0.50  0.5];     % Oil
    
F_boil_r=F_boil./repmat(sum(F_boil,1),size(F_boil,1),1);

P_chp2=P_chp-P_chp1;
F_chp2_r=P_chp2./repmat(sum(P_chp2,1),size(P_chp2,1),1);
%------------------------------------------------------------------------------------------
% 4) Costs and taxes

Cost_r=[0.2 0.2 0.2 0.2 0.2 0.2];              % Cost ranges for constructing supply curve, for example 0.2 means the cost range will be from the "base" until 20% higher. 20% (leave them zero if not needed)
Cost_seg=[5 5 5 5 5 5];                    % Number of cost segments (in constructing supply curve), for example 5 means that the supply curve for that cost range will be divided by 5 (leave them zero if not needed)

carbon_p=8;                                    % Carbon price €/tonne CO2
de_price= 14;                                   % The price in Germany's area

DH_p=[55 60 62 64 60 60];                         % Price of heat from CHP in each area (paid to the plant owner) €/MWh_th

%---------------------------------------------------------------------------------------------------------------------
% 4.1) Marginal cost of power-only plants in each region ("base" cost)

% 4.1.1) Fuel costs (€/MWh_fuel) (for power production)

Cost_fuel=     [0    0    0    0    0   0                % 1. solar PV
                0    0    0    0    0   0                % 2. Onshore wind
                0    0    0    0    0    0               % 3. Offshore wind
                0    0    0    0    0    0              % 4. River hydro (unregulated)
                1    1    1    1    1    1               % 5. Hydro (IMPORTANT: the cost will be later modified based on water value)
                3    3    3    3    3    3               % 6. Nuclear
                0    0    0    0    0   0                % 7. Waste
                14   14   14   14   14   14               % 8. Peat
                22   22   22   22   22    22             % 9. Biomass thermal
                12   12   12   12   12    12             % 10. Coal thermal
                33   33   33   33   33    33             % 11. Gas thermal
                36   36   36   36   36    36             % 12. Oil thermal
                36   36   36   36   36   36];            % 13. Reserve plants (for emergency cases)

% 4.1.2) Variable O&M costs (€/MWh_e)  (only power production), also used for only-power from CHP
Cost_vom=      [0.1    0.1    0.1    0.1    0.1    0.1               % 1. solar PV
                0.5    0.5    0.5    0.5    0.5    0.5               % 2. Onshore wind
                0.6    0.6    0.6    0.6    0.6    0.6               % 3. Offshore wind
                0.1    0.1    0.1    0.1    0.1    0.1               % 4. River hydro (unregulated)
                1      1      1      1      1       1              % 5. Hydro
                0.3    0.3    0.3    0.3    0.3     0.3              % 6. Nuclear
                6      6      6      6      6       6              % 7. Waste
                2      2      2      2      2      2              % 8. Peat
                2      2      2      2      2      2              % 9. Biomass thermal
                0.5    0.5    0.5    0.5    0.5     0.5             % 10. Coal thermal
                0.5    0.5    0.5    0.5    0.5     0.5             % 11. Gas thermal
                0.7    0.7    0.7    0.7    0.7      0.7            % 12. Oil thermal
                1      1      1      1      1  1];                 % 13. Reserve plants (for emergency cases)            
            
% 4.1.3) Additional costs related to energy taxes and other costs not reflected above (€/MWh_fuel)  (only power production)
Cost_add=      [0    0    0    0    0    0              % 1. solar PV
                0    0    0    0    0    0               % 2. Onshore wind
                0    0    0    0    0    0               % 3. Offshore wind
                0    0    0    0    0     0              % 4. River hydro (unregulated)
                0    0    0    0    0     0              % 5. Hydro
                0    0    0    0    0     0              % 6. Nuclear
                0    0    0    0    0     0              % 7. Waste
                0    0    0    0    0    0               % 8. Peat
                0    0    0    0    0     0              % 9. Biomass thermal
                0    0    0    0    0    0               % 10. Coal thermal
                0    0    0    0    0    0               % 11. Gas thermal
                0    0    0    0    0    0               % 12. Oil thermal
                200  200  200  200  200  200];               % 13. Reserve plants (for emergency cases, the costs are reflecting low capacity usage of the plants, therefore lon-term marginal costs)                                  
%--------------------------------------------------------
% 4.1.4) Hydro water value
Wv_seg=[0.3 0.2 0.2 0.2 0.1];               % NO RELATION WITH NUMBER OF regions. Water value main segments
Cost_wv=[1 1.1 1.25 1.5 1.75];              % NO RELATION WITH NUMBER OF regions. The cost increment attributed to the water value segments above

Cost_hyd_coff= [2.85   1.85    0  0   2.    2  % The storage-level dependant cost coefficients for each country (for 7 storage level)
               2.4     1.75    0  0   1.75  2
               1.45    1.45    0  0   1.5   2
               1.25    1.25    0  0   1.25  2
               1.00    1.0     0  0   1.0   2
               0.6     0.75    0  0   0.75  2
               0.5     0.6     0  0   0.6   2];     
     
Cost_hyd_corr=f_watervalue;
%------------------------------------------------------------------------------------------------------------------
% 4.2) Marginal cost of CHP plants

% 4.2.1) Cost of fuels for heat production in CHP (€/MWh_fuel)(including taxes and other additions) 
Cost_fuel_chp1=     [0    0    0    0    0                   % 1. solar PV
                    0    0    0    0    0                   % 2. Onshore wind
                    0    0    0    0    0                   % 3. Offshore wind
                    0    0    0    0    0                   % 4. River hydro (unregulated)
                    0    0    0    0    0                   % 5. Hydro
                    0    0    0    0    0                   % 6. Nuclear
                    0    0    0    0    0                   % 7. Waste
                    17   17   17   17   17                 % 8. Peat
                    22   22   22   22   22                 % 9. Biomass thermal
                    23   23   23   23   23                 % 10. Coal thermal
                    41   41   41   41   41                 % 11. Gas thermal
                    41   41   41   41   41                 % 12. Oil thermal
                    0    0    0    0    0    ];              % 13. Reserve plants (for emergency cases)
Cost_fuel_chp1=[Cost_fuel_chp1 Cost_fuel_chp1(:,5)];   
                
Cost_fuel_chp2=     [0    0    0    0    0                   % 1. solar PV
                    0    0    0    0    0                   % 2. Onshore wind
                    0    0    0    0    0                   % 3. Offshore wind
                    0    0    0    0    0                   % 4. River hydro (unregulated)
                    0    0    0    0    0                   % 5. Hydro
                    0    0    0    0    0                   % 6. Nuclear
                    0    0    0    0    0                   % 7. Waste
                    17   17   17   17   17                 % 8. Peat
                    22   22   22   22   22                 % 9. Biomass thermal
                    23   23   23   23   23                 % 10. Coal thermal
                    41   41   41   41   41                 % 11. Gas thermal
                    41   41   41   41   41                 % 12. Oil thermal
                    0    0    0    0    0    ];              % 13. Reserve plants (for emergency cases)
Cost_fuel_chp2=[Cost_fuel_chp2 Cost_fuel_chp2(:,5)];               
% 4.2.2) Variable O&M costs of CHP plants(€/MWh_e)  (total costs based on power output)
Cost_vom_chp=      [0    0    0    0    0                   % 1. solar PV
                    0    0    0    0    0                    % 2. Onshore wind
                    0    0    0    0    0                    % 3. Offshore wind
                    0    0    0    0    0                    % 4. River hydro (unregulated)
                    0    0    0    0    0                    % 5. Hydro 
                    0    0    0    0    0                     % 6. Nuclear
                    10   10   10   10   10                   % 7. Waste
                    3    3    3    3    3                    % 8. Peat
                    3    3    3    3    3                     % 9. Biomass thermal
                    0.9   0.9   0.9   0.9   0.9             % 10. Coal thermal
                    0.9   0.9   0.9   0.9   0.9             % 11. Gas thermal
                    1    1    1    1    1                     % 12. Oil thermal
                    0    0    0    0    0   ];                % 13. Reserve plants (for emergency cases)
 Cost_vom_chp=[Cost_vom_chp Cost_vom_chp(:,5)];   
%______________________________________________________________________________________________________________________
% 4.3) Market couplings
Market_exch=table2array(A.Market_exch);      % external exchange with the order {'FI_RU','FI_EE','SE_DE','DK1_DE','DK2_DE','NO_NL'}
Exp_uk_2014=load('uk_export_2014.txt');
Market_exch=[Market_exch Exp_uk_2014];
% 5) Pre-calculations

n_chp=size(P_chp,1);                        % Number of CHP production modes

% 5.1) Marginal cost of production plants
% 5.1.1) Marginal cost of power-only plants

Cost_p=round(Cost_fuel./Eff_el+Cost_vom+Cost_add./Eff_el+(Emm_fac./Eff_el).*carbon_p,1);                   % €/MWh_e

% 5.1.2) Marginal cost of CHP plants
% Marginal cost ONLY related to power production from CHP 
Cost_el_chp1=Cost_fuel./(Eff_el_chp1+Eff_th_chp1)+Cost_vom_chp+(Emm_fac./(Eff_el_chp1+Eff_th_chp1)).*carbon_p;    % Local CHP €/MWh_e
Cost_el_chp2=Cost_fuel./(Eff_el_chp2+Eff_th_chp2)+Cost_vom_chp+(Emm_fac./(Eff_el_chp2+Eff_th_chp2)).*carbon_p;    % Central CHP€/MWh_e

% Marginal cost ONLY related to heat production from CHP as a result of ONE unit of power output 
Cost_th_chp1=Cost_fuel_chp1.*(Eff_th_chp1./Eff_el_chp1)./(Eff_el_chp1+Eff_th_chp1)+...
    Emm_fac.*carbon_p.*(Eff_th_chp1./Eff_el_chp1)./(Eff_el_chp1+Eff_th_chp1);        % €/MWh_e

Cost_th_chp2=Cost_fuel_chp2.*(Eff_th_chp2./Eff_el_chp2)./(Eff_el_chp2+Eff_th_chp2)+...
    Emm_fac.*carbon_p.*(Eff_th_chp2./Eff_el_chp2)./(Eff_el_chp2+Eff_th_chp2);        % €/MWh_e

% TOTAL marginal costs of CHP after producing one unit of electriicty 
Cost_chp1_tot=Cost_el_chp1+Cost_th_chp1;                           % €/MWh_e
Cost_chp2_tot=Cost_el_chp2+Cost_th_chp2;                           % €/MWh_e

% Final marginal costs of CHP after deducting revenues from heat sales for CHP
Cost_chp1_mc=Cost_chp1_tot-repmat(DH_p,size(Cost_fuel_chp1,1),1).*(Eff_th_chp1./Eff_el_chp1);           % Marginal cost of CHP plants, potential price offer to the power market (€/MWh_e)
Cost_chp2_mc=Cost_chp2_tot-repmat(DH_p,size(Cost_fuel_chp2,1),1).*(Eff_th_chp2./Eff_el_chp2);

%----------------------------------------------------------------------------------------------------------------------------
% 5.2) Caculation of CHP heat and power cpacity used, and DH demand met by CHP

% Dh_chp_max=Dh.*Dh_chp_r;                                              % Maximum annual DH demand to be met by CHP in each area (MWh_th/h)
Dd_chp1_max=Dd_chpdh.*Dd_chp1_r;                                        % Maximum DH demand to be met by local CHP in each area (MWh_th/h)
Dd_chp2_max=Dd_chpdh-Dd_chp1_max;                                       % Maximum DH demand to be met by central CHP in each area (MWh_th/h)

Dd_chp1_hh=Dd_chp_hh.*repmat(Dd_chp1_r,nhp,1);                           % Hourly DH demand for small CHP plants (MWh_th/h)
Dd_chp2_hh=Dd_chp_hh-Dd_chp1_hh;


T_chp1=P_chp1.*Eff_th_chp1(7:n_chp+6,:)./Eff_el_chp1(7:n_chp+6,:);              % Maximum heat capacity of local CHP plants  (MW_th)
F_chp1_e=P_chp1./repmat(sum(P_chp1,1),n_chp,1);                                 % Share of fuels in electricity production local CHP plants in each area
F_chp1_t=T_chp1./repmat(sum(T_chp1,1),n_chp,1);


T_chp2=P_chp2.*Eff_th_chp2(7:n_chp+6,:)./Eff_el_chp2(7:n_chp+6,:);              % Maximum heat capacity of central CHP plants  (MW_th)
F_chp2_e=P_chp2./repmat(sum(P_chp2,1),n_chp,1);                                 % Share of fuels in electricity production central CHP plants in each area
F_chp2_t=T_chp2./repmat(sum(T_chp2,1),n_chp,1);

Cost_ppchp2_mc=round(Cost_fuel./Eff_el+Cost_vom_chp+Cost_add./Eff_el+(Emm_fac./Eff_el).*carbon_p,1);            % ONLY based on power, assuming the same efficiency as a power plant                  % €/MWh_e

%--------------------------------------------------------------------------------------------------------------------------
% 6) Final results based on number of regions

Name_r=Name_r(1:n_tot);

Eff_el_chp1=Eff_el_chp1(7:size(P_chp,1)+6,:);
Eff_th_chp1=Eff_th_chp1(7:size(P_chp,1)+6,:);

Eff_el_chp2=Eff_el_chp2(7:size(P_chp,1)+6,:);
Eff_th_chp2=Eff_th_chp2(7:size(P_chp,1)+6,:);
Ntc=Ntc(1:n_tot,1:n_tot); 

   