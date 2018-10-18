% Energy system simulation module
%--------------------------------------------------------------------------------
classdef f_esm_dd_r
     properties
            Cost_chp1_pool, Cost_chp2_pool, Dh_pro_chp, Dh_pro_boil, Dh_pro_hob,....
            El_chp1_mar, El_chp2_mar, El_ppchp2_mar, El_ppchp_mar, El_dem_tot, El_dem_res_tot,...
            El_pro_chp1, El_pro_chp, El_pro_nuc, El_pro_vre1, El_pro_vre2, El_pro_vre3, El_pro_vre4,...
            El_pro_vre, Ind_hyd, Name_pool
     end
     methods
         function obj = f_esm_dd_r(n_sys,th_1,th_n)
            
            % 1) Input parameters and loading input data class
            nh=th_n-th_1+1;        % Number of hours in the analysis
            id = f_inputdata_r('Input_data.mat',n_sys);
            %__________________________________________________________________________
            % 2) Initial modifications and preprations
            % CHP (power) capacity used in industry
            P_chp2=id.P_chp-id.P_chp1;

            P_ind_chp1_max=(sum(id.P_chp1)./sum(id.P_chp)).*mean(id.El_chp_ind_hh);
            P_ind_chp2_max=(sum(P_chp2)./sum(id.P_chp)).*mean(id.El_chp_ind_hh);

            nh_p=size(id.El_dem_hh,1);        % Number of hours in input data files

            Dh_dem_tot=id.Dd_hob_hh+id.Dd_chp_ind_hh+id.Dd_chp1_hh+id.Dd_chp2_hh;
            r_chp12=sum(id.P_chp1)./sum(id.P_chp);

            T_chp1_boil=r_chp12.*id.T_chp_boil;
            T_chp2_boil=id.T_chp_boil-T_chp1_boil;

            El_rivhyd_r=id.El_rivhyd_pro./id.P_pow(4,5);
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
            %__________________________________________________________________________
            for r=1:n_sys

            % Heat demand
            dh_sol_hob=id.Dp_sol_hob(r);                         % TWh/a
            dh_sol_chp1=id.Dp_sol_chp1(r);
            dh_sol_chp2=id.Dp_sol_chp2(r);

            st_sol_hob=id.St_sol_hob(r);                       % MWh
            st_sol_chp1=id.St_sol_chp1(r);
            st_sol_chp2=id.St_sol_chp2(r);

            r_sol_hob=id.Dp_sol_hob_r(r);
            r_sol_chp1=id.Dp_sol_chp1_r(r);
            r_sol_chp2=id.Dp_sol_chp2_r(r);

            loss_sol_hob=id.Loss_sol_hob(r);
            loss_sol_chp1=id.Loss_sol_chp1(r);
            loss_sol_chp2=id.Loss_sol_chp2(r);

            Res_sol_th=id.Pv_hh(:,r);   % Ratio

            % DH and CHP
            p_chp1=sum(id.P_chp1(:,r))-P_ind_chp1_max(r);

            eff_el_chp1=mean(id.Eff_el_chp1(:,r),1);
            eff_th_chp1=eff_el_chp1./id.Pth_chp(r);
            p_hp_chp1=id.P_hp_chp1(r);
            cop_hp_chp1=id.Cop_hp_chp1(r);
            t_boil_chp1=T_chp1_boil(r);
            eff_boil_chp1=id.Eff_chp_boil(r);
            st_chp1=id.St_chp1(r);
            r_fixboil_chp1=id.Chp_boil_fix(r);

            p_chp2=sum(P_chp2(:,r))-P_ind_chp2_max(r);
            eff_el_chp2=mean(id.Eff_el_chp2(:,r),1);
            eff_th_chp2=eff_el_chp2./id.Pth_chp(r);
            p_hp_chp2=id.P_hp_chp2(r);
            cop_hp_chp2=id.Cop_hp_chp2(r);
            t_boil_chp2=T_chp2_boil(r);
            eff_boil_chp2=id.Eff_chp_boil(r);
            st_chp2=id.St_chp2(r);
            r_fixboil_chp2=id.Chp_boil_fix(r);

            eff_hob=id.Eff_hob(r);
            %__________________________________________________________________________
            % 3) ANALYSIS
            % 3.1) Heat demand 
            Dh_dem_hob=id.Dd_hob_hh(th_1:th_n,r);
            Dh_dem_chp1=id.Dd_chp1_hh(th_1:th_n,r);
            Dh_dem_chp2=id.Dd_chp2_hh(th_1:th_n,r);

            % Group 1------------------------------------------------------------------
            Stor_sol_hob=zeros(nh,1);
            Dh_p_sol_hob=Res_sol_th(th_1:th_n).*((dh_sol_hob)*1e6/sum(Res_sol_th));

            % Initialization group 1
            if Dh_dem_hob(1)>Dh_p_sol_hob(1)*r_sol_hob+0
                Dh_pro_hob(1,r)=(Dh_dem_hob(1)-Dh_p_sol_hob(1)*r_sol_hob);
            else
                Stor_sol_hob(1)=min(Dh_p_sol_hob(1)*r_sol_hob-Dh_dem_hob(1)+0,st_sol_hob);
            end

            Dh_dem_hob(1)=0;
            for t=2:nh
                if Dh_dem_hob(t)>Dh_p_sol_hob(t)*r_sol_hob+Stor_sol_hob(t-1)*loss_sol_hob
                    Dh_pro_hob(t,r)=Dh_dem_hob(t)-(Dh_p_sol_hob(t)*r_sol_hob+...
                        Stor_sol_hob(t-1)*loss_sol_hob);
                else
                    Stor_sol_hob(t)=min(Dh_p_sol_hob(t)*r_sol_hob+Stor_sol_hob(t-1)*...
                        loss_sol_hob-Dh_dem_hob(t),st_sol_hob);
                end
                 Dh_dem_hob(t)=0;
            end

            % Group 2------------------------------------------------------------------
            Dh_p_sol_chp1=Res_sol_th(th_1:th_n).*((dh_sol_chp1)*1e6/sum(Res_sol_th));

            % Initialization group 2 for solar
            Dh_pro_boil_chp1(:,r)=r_fixboil_chp1.*Dh_dem_chp1;
            Dh_dem_chp1=Dh_dem_chp1-Dh_pro_boil_chp1(:,r);

            % For solar group 2
            if Dh_dem_chp1(1)>Dh_p_sol_chp1(1)*r_sol_chp1+0
                Dh_dem_chp1(1)=Dh_dem_chp1(1)-Dh_p_sol_chp1(1)*r_sol_chp1;
            else
                Stor_sol_chp1(1,r)=min(Dh_p_sol_chp1(1)*r_sol_chp1-Dh_dem_chp1(t)+0,st_sol_chp1);
                Dh_dem_chp1(1)=0;
            end
            for t=2:nh
                if Dh_dem_chp1(t)>Dh_p_sol_chp1(t)*r_sol_chp1+Stor_sol_chp1(t-1,r)*loss_sol_chp1
                    Dh_dem_chp1(t)=Dh_dem_chp1(t)-(Dh_p_sol_chp1(t)*r_sol_chp1+...
                        Stor_sol_chp1(t-1,r)*loss_sol_chp1);
                else
                    Stor_sol_chp1(t,r)=min((Dh_p_sol_chp1(t)*r_sol_chp1+...
                        Stor_sol_chp1(t-1,r)*loss_sol_chp1)-Dh_dem_chp1(t),st_sol_chp1);
                    Dh_dem_chp1(t)=0;
                end
            end

            % Small CHP---------------------------------------------------------------- 
            if eff_el_chp1>0
                for t=1:nh         
                     if Dh_dem_chp1(t)>p_chp1*(eff_th_chp1/eff_el_chp1)                      % 1. CHP
                        Dh_dem_chp1(t)=Dh_dem_chp1(t)-p_chp1*(eff_th_chp1/eff_el_chp1);
                        Dh_pro_chp1(t,r)=p_chp1*(eff_th_chp1/eff_el_chp1);
                        if Dh_dem_chp1(t)>p_hp_chp1*cop_hp_chp1                              % 2. HP (here HP comes after CHP, however 
                             Dh_dem_chp1(t)=Dh_dem_chp1(t)-p_hp_chp1*cop_hp_chp1;               
                             Dh_pro_hp_chp1(t,r)=p_hp_chp1*cop_hp_chp1;                                                                      
                        else
                             Dh_pro_hp_chp1(t,r)=Dh_dem_chp1(t);
                             Dh_dem_chp1(t)=0;
                        end
                        if Dh_dem_chp1(t)>t_boil_chp1-Dh_pro_boil_chp1(t,r)                  % 3. Boiler
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

            % Large CHP----------------------------------------------------------------
            Dh_p_sol_chp2=Res_sol_th(th_1:th_n).*((dh_sol_chp2)*1e6/sum(Res_sol_th));

            % Initialization large CHP
            Dh_pro_boil_chp2(:,r)=r_fixboil_chp2.*Dh_dem_chp2;
            Dh_dem_chp2=Dh_dem_chp2-Dh_pro_boil_chp2(:,r);

            % For solar DH 
            if Dh_dem_chp2(1)>Dh_p_sol_chp2(1)*r_sol_chp2+0
                Dh_dem_chp2(1)=Dh_dem_chp2(1)-Dh_p_sol_chp2(1)*r_sol_chp2;
            else
                Stor_sol_chp2(1,r)=min(Dh_p_sol_chp2(1)*r_sol_chp2-Dh_dem_chp2(t)+0,st_sol_chp2);
                Dh_dem_chp2(1)=0;
            end
            for t=2:nh
                if Dh_dem_chp2(t)>Dh_p_sol_chp2(t)*r_sol_chp2+Stor_sol_chp2(t-1,r)*loss_sol_chp2
                    Dh_dem_chp2(t)=Dh_dem_chp2(t)-(Dh_p_sol_chp2(t)*r_sol_chp2+...
                        Stor_sol_chp2(t-1,r)*loss_sol_chp2);
                else
                    Stor_sol_chp2(t,r)=min(Dh_p_sol_chp2(t)*r_sol_chp2+...
                        Stor_sol_chp2(t-1,r)*loss_sol_chp2-Dh_dem_chp2(t),st_sol_chp2);
                    Dh_dem_chp2(t)=0;
                end
            end

            % For large CHP plants----------------------------------------------------------
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
            %--------------------------------------------------------------------------
            dh_pro_hob=sum(Dh_pro_hob)/1e6;    % TWh_th
            dh_pro_chp1=sum(Dh_pro_chp1)/1e6;
            dh_pro_chp2=sum(Dh_pro_chp2)/1e6;
            
            obj.Dh_pro_hob=Dh_pro_hob;
            obj.Dh_pro_chp=Dh_pro_chp1+Dh_pro_chp2;   % Hourly
            Dh_pro_hp=Dh_pro_hp_chp1+Dh_pro_hp_chp2;
            obj.Dh_pro_boil=Dh_pro_boil_chp1+Dh_pro_boil_chp2;
            Dh_pro_boil_tot=Dh_pro_hob+obj.Dh_pro_boil;
            Dh_pro_sol=Dh_pro_sol_hob+Dh_pro_sol_chp1+Dh_pro_sol_chp2;

            Dh_dem_rem=Dh_dem_hob_rem+Dh_dem_chp1_rem+Dh_dem_chp2_rem;       % The remaining DH demand unfulfilled with the production plants
            %--------------------------------------------------------------------------
            % 3.2) Electricity demand and production
            El_dem_hp_chp1=Dh_pro_hp_chp1./repmat(id.Cop_hp_chp1,nh,1);      % Corresponding electricity demand of heat pumps (in addition to input data of electricity demand)
            El_dem_hp_chp2=Dh_pro_hp_chp2./repmat(id.Cop_hp_chp2,nh,1);

            El_dem_hp_chp=El_dem_hp_chp1+El_dem_hp_chp2;
            El_dem_hp_tot=El_dem_hp_chp;
            obj.El_dem_tot=id.El_dem_hh(th_1:th_n,:)+El_dem_hp_tot;

            % Production---------------------------------------------------------------                        
            obj.El_pro_nuc=id.Nuc_hh(th_1:th_n,:).*repmat(id.P_pow(6,:),nh,1);          % Nuclear hourly production
            el_pro_res1a=id.P_pow(1,:).*id.Cf_pv;                                   % Annual values (MWh/a), initial estimation
            el_pro_res2a=id.P_pow(2,:).*id.Cf_wind*nh_p;
            el_pro_res3a=id.P_pow(3,:).*id.Cf_windoff*nh_p;

            obj.El_pro_vre1=id.Pv_hh(th_1:th_n,:).*repmat(el_pro_res1a./sum(id.Pv_hh,1),nh,1);   % Hourly production MWh/h
            obj.El_pro_vre2=(repmat(id.P_pow(2,:),nh,1).*id.Wind_hh(th_1:th_n,:)).*...
                repmat(id.P_pow(2,:)*nh_p.*id.Cf_wind./(el_pro_res2a+1),nh,1);
            obj.El_pro_vre3=(repmat(id.P_pow(3,:),nh,1).*id.Wind_hh(th_1:th_n,:)).*...
                repmat(id.P_pow(3,:)*nh_p.*id.Cf_windoff./(el_pro_res3a+1),nh,1);
            obj.El_pro_vre4=Hyd_r(th_1:th_n,:).*repmat(id.P_pow(4,:),nh,1);

            obj.El_pro_vre=obj.El_pro_vre1+obj.El_pro_vre2+obj.El_pro_vre3+obj.El_pro_vre4;

            %--------------------------------------------------------------------------
            % 3.3) Power dispatch (technical)
            obj.El_pro_chp1=Dh_pro_chp1.*repmat(id.Pth_chp,nh,1);           % MWh/h Power from CHP based on the respective heat demand
            El_pro_chp2=Dh_pro_chp2.*repmat(id.Pth_chp,nh,1);
            obj.El_pro_chp=obj.El_pro_chp1+El_pro_chp2;

            El_cap_ppchp2=repmat(sum(P_chp2,1),nh,1)-El_pro_chp2;        % Remaining CHP capacity sent to the power market (power only)

            El_dem_res_fix=obj.El_dem_tot-id.El_pro_fix_hh(th_1:th_n,:);     % Residual demand after inflexible production, before RES
            El_dem_res_vre=El_dem_res_fix-obj.El_pro_vre;                    % After RES, excluding hydro storage (and after all inflexible now)
            El_dem_res_nuc=El_dem_res_vre-obj.El_pro_nuc;

            % Later we decide which CHP can be adjusted if a lot of VRE
            El_dem_res_chp1=El_dem_res_nuc-obj.El_pro_chp1;                  % Residual demand after power from local CHP
            El_dem_res_chp2=El_dem_res_chp1-El_pro_chp2;                 % Residual demand after power from central CHP

            obj.El_dem_res_tot=El_dem_res_fix;                               % IT CAN BE modified later if extra VRE

            %--------------------------------------------------------------------------------------
            % 4) Constructing power and cost matrix to the pool
            % 4.1) Participation of CHP in the market (NOTICE: based on CHP strategy) 
            obj.Cost_chp1_pool=round(id.Cost_chp1_mc(7:size(id.P_chp,1)+6,:),1);
            obj.Cost_chp2_pool=round(id.Cost_chp2_mc(7:size(id.P_chp,1)+6,:),1);

            El_chp1_mar=zeros(size(id.P_chp,1),n_sys,nh);
            El_chp2_mar=zeros(size(id.P_chp,1),n_sys,nh);

            El_ppchp2_mar=zeros(size(id.P_chp,1),n_sys,nh);
            El_pp_mar=zeros(size(id.P_pow(7:end,:),1),n_sys,nh);

            for i=1:size(id.P_chp,1)
                El_chp1_mar(i,:,:)=(obj.El_pro_chp1.*repmat(id.P_chp1(i,:)./sum(id.P_chp1),nh,1))';
                El_chp2_mar(i,:,:)=(El_pro_chp2.*repmat(P_chp2(i,:)./sum(P_chp2),nh,1))';

                El_ppchp2_mar(i,:,:)=(El_cap_ppchp2.*repmat(P_chp2(i,:)./sum(P_chp2),nh,1))';
                El_pp_mar(i,:,:)=repmat(id.P_pow(i+6,:),nh,1)';
            end
                El_pp_mar(7,:,:)=repmat(id.P_pow(13,:),nh,1)';

            El_ppchp_mar=zeros(size(El_pp_mar));
            obj.El_ppchp_mar(1:6,:,:)=El_pp_mar(1:6,:,:)+El_ppchp2_mar;       % The sum of power plants and condensing part of unused CHP2
            obj.El_ppchp_mar(7,:,:)= El_pp_mar(7,:,:);   
            obj.El_ppchp2_mar=El_ppchp2_mar;
            
            obj.El_chp1_mar=El_chp1_mar+repmat(id.Cost_chp12_r,size(id.P_chp,1),1,nh).*El_chp2_mar;
            obj.El_chp2_mar=El_chp2_mar-repmat(id.Cost_chp12_r,size(id.P_chp,1),1,nh).*El_chp2_mar;

            obj.Ind_hyd=5:5+length(id.Wv_seg)-1;     % Important: the position of hydro cost segment in the matrix of costs sent to the market            
            obj.Name_pool=cell(size(id.P_pow,1)-1+size(id.P_chp1,1)+size(P_chp2,1)+length(obj.Ind_hyd),1);
            Name= [id.Name_p(1:6),...
                         'CHP1 waste', 'CHP1 peat','CHP1 bio', 'CHP1 coal','CHP1 gas', 'CHP1 oil',...
                         'CHP2 waste', 'CHP2 peat','CHP2 bio', 'CHP2 coal','CHP2 gas', 'CHP2 oil',...
                         id.Name_p(7:end)];
            obj.Name_pool(1:4)=Name(1:4);
            for i=obj.Ind_hyd
                obj.Name_pool{i}=['Hydro' num2str(i-4)];
            end
            obj.Name_pool(5+length(obj.Ind_hyd):end)=Name(6:end);
        end
    end
end
