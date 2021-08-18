% ave_b = HD_exponential(CK,N,S,A,PO,sf)
% * CK:     Temperature Unit: 1 -> Celcius, 0 -> Kelvin
% * N:      Number of strain data points to be used for segmenting the
%           stress-strain curves
% * S:      Specified strain for collecting data for, e.g., plotting
% * A:      Averaging flag for parameter n. 1 -> yes, 0 -> no
% * PO:     Polynomial Order: acceptable values are 'poly1' to 'poly9'.
% * sf:     saving flag. 1 -> yes, 0 -> no
%
% For example: HD_exponential(1,100,0.05,0,'poly9',1)
% and the example show the default values for the function's arguments.
% 
% HD_exponential reads a hot deformation dataset and fit an exponential
% equation on it, and saves all calculations in
% 'MAT_HD_exponential_fit.mat'.
%
% strate*exp(Q/RT) = A*exp(b*sig)
%
% The output includes calculated parameters and their approximations by
% polynomials. The stress values are also available. Comparisons can be
% made for:
% (measured)<-> (predicted)
% ave_b     <-> poly_b
% ave_Q     <-> poly_Q
% ave_lnA   <-> poly_lnA
% stress_m  <-> stress_p
%
% FORMAT OF DATA SET:
% the data files should be stored as "*.data" and their names have to be 
% numbers, e.g., "1.data", "2.data", etc. They are simple text files with 
% the following format: ("LINE #:" should be omitted.)
% LINE 1: temperature (K)
% LINE 2: strain-rate (1/s)
% LINE 3: strain      stress (MPa)
% The format of the rest must be the same as LINE 3.
%
% Dependence: this function needs the following ones to work:
% -HD_g_load_data_files.m
% -HD_g_sortrows_dist.m
% -HD_g_error_cal.m
% 
%% DISCLAIMER
% This program is provided as is for free use.
%
%           Soheil Solhjoo
%           February 17, 2021
% -------------------------------------------------------------------------
function ave_b = HD_exponential(varargin)
%# Check the input argument
idx = ~cellfun('isempty',varargin);
Defaults = {1,100,0.05,0,'poly9',1};
Defaults(idx) = varargin(idx);
[CK,N,strain_spec,averaging,p_order,saving_f] = Defaults{:};
R = 8.314;  %Gas constant
%% Load *.data files
[ss_approx,strain,strate,temp,u_rates,u_temp,index] = ...
    HD_g_load_data_files(N,strain_spec,CK,0);
%% Fitting in steps: strate*exp(Q/RT) = A exp(b.sig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over temp for estimating n
% For that, the stress curves should be sorted for deformation temperature.
[temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx'],1);
% The loop over temp
[~,~,temp(:,2)] = unique(temp);
l_temp = 1;
j = 1;
while l_temp
    cond_temp = temp(:,2) == temp(l_temp,2);
    t_cond_temp = cond_temp;
    t_cond_temp(t_cond_temp==0) = [];
    num = size(t_cond_temp,1);
    
    t_rate = strate(l_temp:l_temp+num-1,1);
    t_stress = ss_approx(l_temp:l_temp+num-1,:);
    
    %# Collect data for figure at the selected strain
    for i=index
        t_X = t_stress(:,i);
        curves_b(:,j) = [temp(l_temp,1) polyfit(t_X,log(t_rate),1) ...
            t_X' log(t_rate')];
    end
    
    %# Measure all data
    for i=1:N
        t_X = t_stress(:,i);
        p_b(:,i) = polyfit(t_X,log(t_rate),1);
    end
    b(j,:) = p_b(1,:);
    
    l_temp = l_temp + num;
    j = j+1;
    if cond_temp(end)==1
        clear p_n
        l_temp = 0;
    end
    clear p_b
end
temp = temp(:,1);
ave_b = mean(b);
c_b   = flip(coeffvalues(fit(strain,ave_b',p_order)));
if averaging
    ave_b = mean(b(:)) * ones(1,size(b,2));
end
clear t_temp l_temp t_rate t_stress num i j cond_temp t_cond_temp t_X p_var
%# Loop over strate for estimating Q
%# The stress curves should be sorted for deformation strain rate.
[temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx],2);
%# The loop over strate
[~,~,strate(:,2)] = unique(strate);
l_strate = 1;
j = 1;
while l_strate
    cond_strate = strate(:,2) == strate(l_strate,2);
    t_cond_strate = cond_strate;
    t_cond_strate(t_cond_strate==0) = [];
    num = size(t_cond_strate,1);
    
    t_temp = temp(l_strate:l_strate+num-1,1)';
    t_stress = ss_approx(l_strate:l_strate+num-1,:);
    
    %# Collect data for figure at the selected strain
    for i=index
        t_X = 1000./(R*t_temp);
        t_Y = ave_b(:,i).*t_stress(:,i);
        curves_Q(:,j) = [strate(l_strate,1) polyfit(t_X,t_Y',1) ...
            t_X t_Y'];
    end
    
    %# Measure all data
    for i=1:N
        t_X = 1000./(R*t_temp);
        t_Y = ave_b(:,i).*t_stress(:,i);
        p_Q(:,i) = polyfit(t_X,t_Y',1);
    end
    Q(j,:) = p_Q(1,:);
    lnA(j,:) = log(strate(l_strate,1)) - p_Q(2,:);
    
    l_strate = l_strate + num;
    j = j+1;
    if cond_strate(end)==1
        clear p_Q
        l_strate = 0;
    end
    clear p_Q
end
strate = strate(:,1);
ave_Q = mean(Q)';
ave_lnA = mean(lnA)';

ave_b = ave_b';
%% Fit polynomial to averaged n,Q,lnA
p_ave_b     = fit(strain,ave_b,p_order);
p_ave_Q     = fit(strain,ave_Q,p_order);
p_ave_lnA   = fit(strain,ave_lnA,p_order);

c_Q   = flip(coeffvalues(p_ave_Q));
c_lnA = flip(coeffvalues(p_ave_lnA));

poly_b   = feval(p_ave_b,strain);
poly_Q   = feval(p_ave_Q,strain);
poly_lnA = feval(p_ave_lnA,strain);
%% Stress values from the model
[mesh_strate,mesh_temp,mesh_b,mesh_Q,mesh_lnA,Z] = ...
    HD_g_mesher(strain,strate,temp,poly_b,poly_Q,poly_lnA,R,1000);
sig = log(Z./exp(mesh_lnA))./mesh_b;

Correlation = sortrows([ss_approx(:),sig(:)]);
stress_m = Correlation(:,1);
stress_p = Correlation(:,2);
errors = HD_g_error_cal(stress_m,stress_p);
%% save the workspace
if saving_f
    save('MAT_HD_exponential_fit',...
        'ss_approx','strain','strate','temp','u_rates','u_temp','index',...
        'strain_spec','averaging','p_order',...
        'ave_b','ave_Q','ave_lnA','poly_b','poly_Q','poly_lnA',...
        'c_b','c_Q','c_lnA','curves_b','curves_Q','b','Q','lnA',...
        'mesh_strate','mesh_temp','mesh_Q','mesh_lnA','sig',...
        'stress_m','stress_p','errors');
end
end