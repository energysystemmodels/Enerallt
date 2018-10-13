% Energy system simulation module
%--------------------------------------------------------------------------------
function[Ntc,Market_exch,El_dem_tot,El_dem_res_tot,Eff_el,Emm_fac,Eff_el_chp1,Eff_el_chp2,Eff_chp_boil,Eff_hob,....
    El_pro_vre1,El_pro_vre2,El_pro_vre3,El_pro_vre4,El_pro_nuc,El_chp_ind_hh,El_pro_fix_hh,...
    El_cap_hyd,El_eff_hyd,Hyd_res_dd,Hyd_infl_dd,St_hyd,Cost_hyd_corr,Wv_seg,Cost_wv,Cost_hyd_coff,St_begin,St_min,Ind_hyd,...
    El_pro_chp1,El_pro_chp,Pth_chp,Dh_pro_chp,Dh_pro_boil,Dh_pro_hob,F_chp1_r,F_chp2_r,F_boil_r,F_hob_r,El_chp1_mar,El_chp2_mar,El_ppchp2_mar,El_ppchp_mar,...
    Cost_chp1_pool,Cost_chp2_pool,Cost_p,Cost_r,Cost_seg,Name_pool,Name_r]=f_esm_dd(n_sys,th_1,th_n)

nh=th_n-th_1+1;        % Number of hours in the analysis

[El_dem_hh,P_pow,Eff_el,Emm_fac,El_chp_ind_hh,El_pro_fix_hh,El_pro_nuc_r,...
 Dd_hob_hh,Dd_chp_ind_hh,Dd_chp1_hh,Dd_chp2_hh,Chp_ind_hh,P_chp,P_chp1,Eff_el_chp1,Eff_el_chp2,Eff_th_chp1,Eff_th_chp2,Eff_chp_boil,Eff_hob,...
 Pth_chp,F_hob_r,F_chp1_r,F_chp2_r,F_boil_r,P_hp_chp1,P_hp_chp2,Cop_hp_chp1,Cop_hp_chp2,T_chp_boil,Chp_boil_fix,...
 St_chp1,St_chp2,Dp_sol_hob,Dp_sol_chp1,Dp_sol_chp2,Dp_sol_hob_r,Dp_sol_chp1_r,Dp_sol_chp2_r,St_sol_hob,St_sol_chp1,St_sol_chp2,Loss_sol_hob,Loss_sol_chp1,Loss_sol_chp2,... 
 Cf_pv,Cf_wind,Cf_windoff,Pv_hh,Wind_hh,Hyd_res_hh,Hyd_infl_hh,El_rivhyd_pro,Nuc_hh,Hyd_res_dd,Hyd_infl_dd,St_hyd,Cost_hyd_corr,Wv_seg,Cost_wv,Cost_hyd_coff,St_begin,St_min,...
 Name_r,Name_p,Name_chp,Ntc,Market_exch,Cost_p,Cost_chp1_mc,Cost_chp2_mc,Cost_ppchp2_mc,Cost_r,Cost_seg]=f_inputdata(n_sys);

clear f_id_uk2
%----------------------------------------------------------
P_pow(10,4)=P_pow(10,4)+2000;   % Coal east denmark correction
El_rivhyd_r=El_rivhyd_pro./P_pow(4,5);   
%-------------------------------------------------------------

% Share of CHP2 with cost method as CHP1
Cost_chp12_r=[0.0  0  0  0.0  0 0];

% CHP (power) capacity used in industry
P_chp2=P_chp-P_chp1;

P_ind_chp1_max=(sum(P_chp1)./sum(P_chp)).*mean(El_chp_ind_hh);
P_ind_chp2_max=(sum(P_chp2)./sum(P_chp)).*mean(El_chp_ind_hh);

nh_p=size(El_dem_hh,1);        % Number of hours in input data files

Dh_dem_tot=Dd_hob_hh+Dd_chp_ind_hh+Dd_chp1_hh+Dd_chp2_hh;
r_chp12=sum(P_chp1)./sum(P_chp);

T_chp1_boil=r_chp12.*T_chp_boil;
T_chp2_boil=T_chp_boil-T_chp1_boil;

Hyd_r=repmat(El_rivhyd_r,1,n_sys);    % Ratio MUST BE CORRECTED AFTER RIVER HYDRO FROM FI AND SE
Hyd_r(isnan(Hyd_r))=0;

Dh_pro_hob=zeros(nh,n_sys);                % Hourly
Dh_pro_chp1=zeros(nh,n_sys);
Dh_pro_chp2=zeros(nh,n_sys);

Dh_pro_boil_chp1=zeros(nh,n_sys);
Dh_pro_boil_chp2=zeros(nh,n_sys);

Dh_pro_sol_hob=zeros(nh,n_sys);
Dh_pro_sol_chp1=zeros(nh,n_sys);
Dh_pro_sol_chp2=zeros(nh,n_sys);
Stor_sol_chp1=zeros(nh,n_sys);            % MWh_th/h
Stor_sol_chp2=zeros(nh,n_sys);

Dh_pro_hp_chp1=zeros(nh,n_sys);
Dh_pro_hp_chp2=zeros(nh,n_sys);

Dh_dem_hob_rem=zeros(nh,n_sys);            % Remaining DH demand not met by the fleet!
Dh_dem_chp1_rem=zeros(nh,n_sys);
Dh_dem_chp2_rem=zeros(nh,n_sys);

for r=1:n_sys
    
% Heat demand
% dh_dem_hob=sum(Dd_hob_hh(:,r));
% dh_dem_chp1=sum(Dd_chp1_hh(:,r));
% dh_dem_chp2=sum(Dd_chp2_hh(:,r));

dh_sol_hob=Dp_sol_hob(r);                         % TWh/a
dh_sol_chp1=Dp_sol_chp1(r);
dh_sol_chp2=Dp_sol_chp2(r);

st_sol_hob=St_sol_hob(r);                       % MWh
st_sol_chp1=St_sol_chp1(r);
st_sol_chp2=St_sol_chp2(r);

r_sol_hob=Dp_sol_hob_r(r);
r_sol_chp1=Dp_sol_chp1_r(r);
r_sol_chp2=Dp_sol_chp2_r(r);

loss_sol_hob=Loss_sol_hob(r);
loss_sol_chp1=Loss_sol_chp1(r);
loss_sol_chp2=Loss_sol_chp2(r);
       
Res_sol_th=Pv_hh(:,r);   % Ratio

% DH and CHP
p_chp1=sum(P_chp1(:,r))-P_ind_chp1_max(r);

eff_el_chp1=mean(Eff_el_chp1(:,r),1);
eff_th_chp1=eff_el_chp1./Pth_chp(r);
p_hp_chp1=P_hp_chp1(r);
cop_hp_chp1=Cop_hp_chp1(r);
t_boil_chp1=T_chp1_boil(r);
eff_boil_chp1=Eff_chp_boil(r);
st_chp1=St_chp1(r);
r_fixboil_chp1=Chp_boil_fix(r);

p_chp2=sum(P_chp2(:,r))-P_ind_chp2_max(r);
eff_el_chp2=mean(Eff_el_chp2(:,r),1);
eff_th_chp2=eff_el_chp2./Pth_chp(r);
p_hp_chp2=P_hp_chp2(r);
cop_hp_chp2=Cop_hp_chp2(r);
t_boil_chp2=T_chp2_boil(r);
eff_boil_chp2=Eff_chp_boil(r);
st_chp2=St_chp2(r);
r_fixboil_chp2=Chp_boil_fix(r);

eff_hob=Eff_hob(r);

% B) ANALYSIS

% 1) Heat demand 
Dh_dem_hob=Dd_hob_hh(th_1:th_n,r);
Dh_dem_chp1=Dd_chp1_hh(th_1:th_n,r);
Dh_dem_chp2=Dd_chp2_hh(th_1:th_n,r);

% Group 1
Stor_sol_hob=zeros(nh,1);
% if sum(Res_sol_th(th_1:th_n))>0
    Dh_p_sol_hob=Res_sol_th(th_1:th_n).*((dh_sol_hob)*1e6/sum(Res_sol_th));
% else
%     Dh_pro_sol_g1(th_1:th_n)=0;
% end

% Initialization group 1
if Dh_dem_hob(1)>Dh_p_sol_hob(1)*r_sol_hob+0
    Dh_pro_hob(1,r)=(Dh_dem_hob(1)-Dh_p_sol_hob(1)*r_sol_hob);
else
     Stor_sol_hob(1)=min(Dh_p_sol_hob(1)*r_sol_hob-Dh_dem_hob(1)+0,st_sol_hob);
end
Dh_dem_hob(1)=0;
for t=2:nh
    if Dh_dem_hob(t)>Dh_p_sol_hob(t)*r_sol_hob+Stor_sol_hob(t-1)*loss_sol_hob
     Dh_pro_hob(t,r)=Dh_dem_hob(t)-(Dh_p_sol_hob(t)*r_sol_hob+Stor_sol_hob(t-1)*loss_sol_hob);
    else
     Stor_sol_hob(t)=min(Dh_p_sol_hob(t)*r_sol_hob+Stor_sol_hob(t-1)*loss_sol_hob-Dh_dem_hob(t),st_sol_hob);
    end
     Dh_dem_hob(t)=0;
end

% Group 2

% if sum(Res_sol_th(th_1:th_n))>0
    Dh_p_sol_chp1=Res_sol_th(th_1:th_n).*((dh_sol_chp1)*1e6/sum(Res_sol_th));
% else
%     Dh_pro_sol_g2(th_1:th_n)=0;
% end

% Initialization group 2 for solar
Dh_pro_boil_chp1(:,r)=r_fixboil_chp1.*Dh_dem_chp1;
Dh_dem_chp1=Dh_dem_chp1-Dh_pro_boil_chp1(:,r);

% For solar group 2
if Dh_dem_chp1(1)>Dh_p_sol_chp1(1)*r_sol_chp1+0;
    Dh_dem_chp1(1)=Dh_dem_chp1(1)-Dh_p_sol_chp1(1)*r_sol_chp1;
else
    Stor_sol_chp1(1,r)=min(Dh_p_sol_chp1(1)*r_sol_chp1-Dh_dem_chp1(t)+0,st_sol_chp1);
    Dh_dem_chp1(1)=0;
end
for t=2:nh
    if Dh_dem_chp1(t)>Dh_p_sol_chp1(t)*r_sol_chp1+Stor_sol_chp1(t-1,r)*loss_sol_chp1;
     Dh_dem_chp1(t)=Dh_dem_chp1(t)-(Dh_p_sol_chp1(t)*r_sol_chp1+Stor_sol_chp1(t-1,r)*loss_sol_chp1);
    else
     Stor_sol_chp1(t,r)=min((Dh_p_sol_chp1(t)*r_sol_chp1+Stor_sol_chp1(t-1,r)*loss_sol_chp1)-Dh_dem_chp1(t),st_sol_chp1);
     Dh_dem_chp1(t)=0;
    end
end

% Small CHP 

if eff_el_chp1>0
for t=1:nh         
     if Dh_dem_chp1(t)>p_chp1*(eff_th_chp1/eff_el_chp1)                      % 1. CHP
       Dh_dem_chp1(t)=Dh_dem_chp1(t)-p_chp1*(eff_th_chp1/eff_el_chp1);
       Dh_pro_chp1(t,r)=p_chp1*(eff_th_chp1/eff_el_chp1);
         if Dh_dem_chp1(t)>p_hp_chp1*cop_hp_chp1                                    % 2. HP (here HP comes after CHP, however 
              Dh_dem_chp1(t)=Dh_dem_chp1(t)-p_hp_chp1*cop_hp_chp1;               % I should model the case HP can come before CHP if excess power)
              Dh_pro_hp_chp1(t,r)=p_hp_chp1*cop_hp_chp1;                                                                      
         else
             Dh_pro_hp_chp1(t,r)=Dh_dem_chp1(t);
             Dh_dem_chp1(t)=0;
         end
         if Dh_dem_chp1(t)>t_boil_chp1-Dh_pro_boil_chp1(t,r)              % 3. Boil
             Dh_dem_chp1(t)=Dh_dem_chp1(t)-(t_boil_chp1-Dh_pro_boil_chp1(t,r));
             Dh_pro_boil_chp1(t,r)=t_boil_chp1;
         else
             Dh_pro_boil_chp1(t,r)=Dh_pro_boil_chp1(t,r)+Dh_dem_chp1(t);
             Dh_dem_chp1(t)=0;
         end
     else
      Dh_pro_chp1(t,r)=Dh_dem_chp1(t);
      Dh_dem_chp1(t)=0;
     end
end
end

% Large CHP

Dh_p_sol_chp2=Res_sol_th(th_1:th_n).*((dh_sol_chp2)*1e6/sum(Res_sol_th));

% Initialization large CHP
Dh_pro_boil_chp2(:,r)=r_fixboil_chp2.*Dh_dem_chp2;
Dh_dem_chp2=Dh_dem_chp2-Dh_pro_boil_chp2(:,r);

% For solar DH 
if Dh_dem_chp2(1)>Dh_p_sol_chp2(1)*r_sol_chp2+0;
    Dh_dem_chp2(1)=Dh_dem_chp2(1)-Dh_p_sol_chp2(1)*r_sol_chp2;
else
    Stor_sol_chp2(1,r)=min(Dh_p_sol_chp2(1)*r_sol_chp2-Dh_dem_chp2(t)+0,st_sol_chp2);
    Dh_dem_chp2(1)=0;
end
for t=2:nh
    if Dh_dem_chp2(t)>Dh_p_sol_chp2(t)*r_sol_chp2+Stor_sol_chp2(t-1,r)*loss_sol_chp2;
     Dh_dem_chp2(t)=Dh_dem_chp2(t)-(Dh_p_sol_chp2(t)*r_sol_chp2+Stor_sol_chp2(t-1,r)*loss_sol_chp2);
    else
     Stor_sol_chp2(t,r)=min(Dh_p_sol_chp2(t)*r_sol_chp2+Stor_sol_chp2(t-1,r)*loss_sol_chp2-Dh_dem_chp2(t),st_sol_chp2);
      Dh_dem_chp2(t)=0;
    end
end

% For CHP group 3

if eff_el_chp2>0
for t=1:nh
     if Dh_dem_chp2(t)>p_chp2*(eff_th_chp2/eff_el_chp2)
      Dh_dem_chp2(t)=Dh_dem_chp2(t)-p_chp2*(eff_th_chp2/eff_el_chp2);
      Dh_pro_chp2(t,r)=p_chp2*(eff_th_chp2/eff_el_chp2);
         if Dh_dem_chp2(t)>p_hp_chp2*cop_hp_chp2
             Dh_dem_chp2(t)=Dh_dem_chp2(t)-p_hp_chp2*cop_hp_chp2;
             Dh_pro_hp_chp2(t,r)=p_hp_chp2*cop_hp_chp2;
         else
             Dh_pro_hp_chp2(t,r)=Dh_dem_chp2(t);
             Dh_dem_chp2(t)=0;
         end
         if Dh_dem_chp2(t)>t_boil_chp2-Dh_pro_boil_chp2(t,r)
             Dh_dem_chp2(t)=Dh_dem_chp2(t)-(t_boil_chp2-Dh_pro_boil_chp2(t,r));
             Dh_pro_boil_chp2(t,r)=t_boil_chp2;
         else
             Dh_pro_boil_chp2(t,r)=Dh_pro_boil_chp2(t,r)+Dh_dem_chp2(t);
             Dh_dem_chp2(t)=0;
         end
     else
      Dh_pro_chp2(t,r)=Dh_dem_chp2(t);
      Dh_dem_chp2(t)=0;
     end
end
end
Dh_pro_sol_hob(:,r)=Dh_p_sol_hob;
Dh_pro_sol_chp1(:,r)=Dh_p_sol_chp1;
Dh_pro_sol_chp2(:,r)=Dh_p_sol_chp2;

if find(Dh_dem_hob)
    warning(['DH demand not met by the existing fleet for heat-only boilers in region no.' num2str(r) ' !'])
end
if find(Dh_dem_chp1)
    warning(['DH demand not met by the existing fleet for small CHP in region no.' num2str(r) ' !'])
end
if find(Dh_dem_chp2)
    warning(['DH demand not met by the existing fleet for small CHP in region no.' num2str(r) ' !'])
end
Dh_dem_hob_rem(:,r)=Dh_dem_hob;            % Remaining DH demand not met by the fleet!
Dh_dem_chp1_rem(:,r)=Dh_dem_chp1;
Dh_dem_chp2_rem(:,r)=Dh_dem_chp2;
end
%-------------------------------------------------------------------------------------
dh_pro_hob=sum(Dh_pro_hob)/1e6;    % TWh_th
dh_pro_chp1=sum(Dh_pro_chp1)/1e6;
dh_pro_chp2=sum(Dh_pro_chp2)/1e6;


Dh_pro_chp=Dh_pro_chp1+Dh_pro_chp2;   % Hourly
Dh_pro_hp=Dh_pro_hp_chp1+Dh_pro_hp_chp2;
Dh_pro_boil=Dh_pro_boil_chp1+Dh_pro_boil_chp2;
Dh_pro_boil_tot=Dh_pro_hob+Dh_pro_boil;
Dh_pro_sol=Dh_pro_sol_hob+Dh_pro_sol_chp1+Dh_pro_sol_chp2;

Dh_dem_rem=Dh_dem_hob_rem+Dh_dem_chp1_rem+Dh_dem_chp2_rem;       % The remaining DH demand unfulfilled with the production plants

%---------------------------------------------------------------------------------------
% 2) Electricity demand and production

El_dem_hp_chp1=Dh_pro_hp_chp1./repmat(Cop_hp_chp1,nh,1);    % Corresponding electricity demand of heat pumps (in addition to input data of electricity demand)
El_dem_hp_chp2=Dh_pro_hp_chp2./repmat(Cop_hp_chp2,nh,1);

El_dem_hp_chp=El_dem_hp_chp1+El_dem_hp_chp2;

% el_dem_hp_hh=h_dem_hp_hh/cop_hp_hh;       % Elecricity demand for individual HPs (households)
% El_dem_hp_tot=El_dem_hp_chp+El_dem_hp_hh;
El_dem_hp_tot=El_dem_hp_chp;

% el_dem_a=el_dem+el_dem_ec+el_dem_hp_tot;    % Annual el demand, el_dem_eh is included in el_dem, Add other items later if any
El_dem_tot=El_dem_hh(th_1:th_n,:)+El_dem_hp_tot;

% Production                        
El_pro_nuc=Nuc_hh(th_1:th_n,:).*repmat(P_pow(6,:),nh,1);                     % Nuclear hourly production
el_pro_res1a=P_pow(1,:).*Cf_pv;                           % Annual values (MWh/a), initial estimation
el_pro_res2a=P_pow(2,:).*Cf_wind*nh_p;
el_pro_res3a=P_pow(3,:).*Cf_windoff*nh_p;

El_pro_vre1=Pv_hh(th_1:th_n,:).*repmat(el_pro_res1a./sum(Pv_hh,1),nh,1);   % Hourly production MWh/h
El_pro_vre2=(repmat(P_pow(2,:),nh,1).*Wind_hh(th_1:th_n,:)).*repmat(P_pow(2,:)*nh_p.*Cf_wind./(el_pro_res2a+1),nh,1);
El_pro_vre3=(repmat(P_pow(3,:),nh,1).*Wind_hh(th_1:th_n,:)).*repmat(P_pow(3,:)*nh_p.*Cf_windoff./(el_pro_res3a+1),nh,1);
El_pro_vre4=Hyd_r(th_1:th_n,:).*repmat(P_pow(4,:),nh,1);

El_pro_vre=El_pro_vre1+El_pro_vre2+El_pro_vre3+El_pro_vre4;
%------------------------------------------------------------------------------------------
% 2.2. Power dispatch (technical)

El_pro_chp1=Dh_pro_chp1.*repmat(Pth_chp,nh,1);           % MWh/h Power from CHP based on the respective heat demand
El_pro_chp2=Dh_pro_chp2.*repmat(Pth_chp,nh,1);
El_pro_chp=El_pro_chp1+El_pro_chp2;

El_cap_ppchp2=repmat(sum(P_chp2,1),nh,1)-El_pro_chp2;                  % Remaining CHP capacity sent to the power market (power only)

El_dem_res_fix=El_dem_tot-El_pro_fix_hh(th_1:th_n,:);                  % Residual demand after inflexible production, before RES
El_dem_res_vre=El_dem_res_fix-El_pro_vre;                                % After RES, excluding hydro storage (and after all inflexible now)
El_dem_res_nuc=El_dem_res_vre-El_pro_nuc;

% Later we decide which CHP can be adjusted if a lot of VRE
El_dem_res_chp1=El_dem_res_nuc-El_pro_chp1;                                    % Residual demand after power from local CHP
El_dem_res_chp2=El_dem_res_chp1-El_pro_chp2;                                   % Residual demand after power from central CHP

El_dem_res_tot=El_dem_res_fix;          % IT CAN BE modified later if extra VRE
%--------------------------------------------------------------------------------------
% 6) Constructing power and cost matrix to the pool

% % 6.1) Participation of CHP in the market (NOTICE: based on CHP strategy) 
Cost_chp1_pool=round(Cost_chp1_mc(7:size(P_chp,1)+6,:),1);
Cost_chp2_pool=round(Cost_chp2_mc(7:size(P_chp,1)+6,:),1);
%-------------------------------------------------------------------------------
El_chp1_mar=zeros(size(P_chp,1),n_sys,nh);
El_chp2_mar=zeros(size(P_chp,1),n_sys,nh);

El_ppchp2_mar=zeros(size(P_chp,1),n_sys,nh);
El_pp_mar=zeros(size(P_pow(7:end,:),1),n_sys,nh);

for i=1:size(P_chp,1)
    El_chp1_mar(i,:,:)=(El_pro_chp1.*repmat(P_chp1(i,:)./sum(P_chp1),nh,1))';
    El_chp2_mar(i,:,:)=(El_pro_chp2.*repmat(P_chp2(i,:)./sum(P_chp2),nh,1))';
    
    El_ppchp2_mar(i,:,:)=(El_cap_ppchp2.*repmat(P_chp2(i,:)./sum(P_chp2),nh,1))';
    El_pp_mar(i,:,:)=repmat(P_pow(i+6,:),nh,1)';
end
    El_pp_mar(7,:,:)=repmat(P_pow(13,:),nh,1)';

El_ppchp_mar=zeros(size(El_pp_mar));
El_ppchp_mar(1:6,:,:)=El_pp_mar(1:6,:,:)+El_ppchp2_mar;       % The sum of power plants and condensing part of unused CHP2
El_ppchp_mar(7,:,:)= El_pp_mar(7,:,:);   

El_chp1_mar=El_chp1_mar+repmat(Cost_chp12_r,size(P_chp,1),1,nh).*El_chp2_mar;
El_chp2_mar=El_chp2_mar-repmat(Cost_chp12_r,size(P_chp,1),1,nh).*El_chp2_mar;

El_cap_hyd=P_pow(5,:);
El_eff_hyd=Eff_el(5,:);

Ind_hyd=5:5+length(Wv_seg)-1;                                % Important: the position of hydro cost segment in the matrix of costs sent to the market            
Name_pool=cell(size(P_pow,1)-1+size(P_chp1,1)+size(P_chp2,1)+length(Ind_hyd),1);
Name= [Name_p(1:6),...
             'CHP1 waste', 'CHP1 peat','CHP1 bio', 'CHP1 coal','CHP1 gas', 'CHP1 oil',...
             'CHP2 waste', 'CHP2 peat','CHP2 bio', 'CHP2 coal','CHP2 gas', 'CHP2 oil',...
             Name_p(7:end)];
Name_pool(1:4)=Name(1:4);
for i=Ind_hyd
    Name_pool{i}=['Hydro' num2str(i-4)];
end
Name_pool(5+length(Ind_hyd):end)=Name(6:end);
