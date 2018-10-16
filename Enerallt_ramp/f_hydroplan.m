function [Hyd, Cost_mod, stor_d]=f_hydroplan(Cost_wv,WV_seg,Cost_coff,hyd_reg,stor_pre,st_min,...
                                             hyd_rest,Dem,Vres,el_cap,st_max,eff)

stor_ava=stor_pre+hyd_rest;       % Available hydro in the storage from the previous day
hyd_ava= hyd_reg+stor_ava;    % Available Hydro resources for that day (MWh/d)   

st_min=0.1*hyd_ava;                  % Minimum storage level for each day

if hyd_ava<=0.25*(st_max-st_min) && hyd_ava>=0        % c: Correction factor of prices based on the sales in the previous day
    cost_coff=Cost_coff(1);
elseif   hyd_ava>0.25*(st_max-st_min) && hyd_ava<=0.5*(st_max-st_min)
    cost_coff=Cost_coff(2);
elseif   hyd_ava>0.5*(st_max-st_min)&& hyd_ava<=0.75*(st_max-st_min)
    cost_coff=Cost_coff(3);
elseif   hyd_ava>0.75*(st_max-st_min) && hyd_ava<=(st_max-st_min)
    cost_coff=Cost_coff(4);
elseif   hyd_ava>(st_max-st_min)&& hyd_ava<=st_max
    cost_coff=Cost_coff(5);
elseif   hyd_ava>st_max && hyd_ava<=1.25*st_max
    cost_coff=Cost_coff(6);
    hyd_ava=st_max;                     % The rest of hydro will be spilled
elseif   hyd_ava>1.25*st_max 
    cost_coff=Cost_coff(7);
    hyd_ava=st_max;                     % The rest of hydro will be spilled

elseif fix(hyd_ava)<0
    error('Hydro storage level negative!!')
end

% Water value
Cost_mod=cost_coff.*Cost_wv;    % Coefficients for modification of cost based on the earnings in the prvious day

%-----------------------------------------------------
% Analysis

Residual=Dem-Vres;
Ratio=Residual./sum(Residual);
Ratio(Ratio<0)=0;
Q=zeros(24,1);    % Quantity of hydro sent to the market in each hour (MWh/h)
for i=1:24
    Q(i)=min(Ratio(i)*(hyd_ava)*eff,el_cap);
end
stor_d=hyd_ava-(sum(Q)/eff);  
Hyd=fix(repmat(Q,1,length(WV_seg)).*repmat(WV_seg,24,1));
