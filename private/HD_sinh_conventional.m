% [mesh_alpha,n,Q,lnA] = HD_sinh_conventional(CK,N,S,A,PO,sf):
% * CK:     Temperature Unit: 1 -> Celcius, 0 -> Kelvin
% * N:      Number of strain data points to be used for segmenting the
%           stress-strain curves
% * S:      Specified strain for collecting data for, e.g., plotting
% * A:      Averaging flag for parameters [alpha,n]: 1 -> yes, 0 -> no
% * PO:     Polynomial Order: acceptable values are 'poly1' to 'poly9'.
% * sf:     saving flag. 1 -> yes, 0 -> no
%
% For example: HD_sinh_conventional(1,100,0.05,[1,1],'poly9',1)
% and the example show the default values for the function's arguments.
% 
% HD_sinh_conventional reads a hot deformation dataset and fit Sellars and
% Tegart's hyperbolic sine function on it, and save all calculations in
% 'MAT_HD_sinh_conventional_fit.mat'.
%
% strate*exp(Q/RT) = A*(sinh(alpha*sigma))^n
%
% The output includes calculated parameters and their approximations by
% polynomials. The stress values are also available. Comparisons can be
% made for:
% (measured)<-> (predicted)
% ave_alpha <-> poly_alpha
% ave_n     <-> poly_n
% ave_Q     <-> poly_Q
% ave_lnA   <-> poly_lnA
% stress_m  <-> stress_p
%
% FORMAT OF DATA SET: the data files should be stored as "*.data" and their
% names have to be numbers, e.g., "1.data", "2.data", etc. They are simple
% text files with the following format: ("LINE #:" should be omitted.)
% LINE 1: temperature (K)
% LINE 2: strain-rate (1/s)
% LINE 3: strain      stress (MPa)
% The format of the rest must be the same as LINE 3.
%
% Dependence: this function needs the following ones to work:
% -HD_g_load_data_files.m
% -HD_g_sortrows_dist.m
% -HD_g_error_cal.m
% -HD_power_law.m
% -HD_exponential.m
% 
%% DISCLAIMER
% This program is provided as is for free use.
%
%           Soheil Solhjoo
%           February 17, 2021
% -------------------------------------------------------------------------
function [mesh_alpha,mesh_n,Q,lnA] = HD_sinh_conventional(varargin)
%# Check the input argument
idx = ~cellfun('isempty',varargin);
Defaults = {1,100,0.05,[0,0],'poly9',1};
Defaults(idx) = varargin(idx);
[CK,N,strain_spec,averaging,p_order,saving_f] = Defaults{:};
R = 8.314;  %Gas constant

[averaging_alpha,averaging_n] = feval(@(x) x{:}, num2cell(averaging));
%% Find & Calculate ava_n, ave_b and ave_alpha
ave_n = HD_power_law(CK,100,0,0,'poly9',0);
ave_b = HD_exponential(CK,100,0,0,'poly9',0);

ave_alpha = ave_b./ave_n;

if averaging_alpha
    ave_alpha = ones(size(ave_alpha)) * mean(ave_alpha);
end
%% Load *.data files
[ss_approx,strain,strate,temp,u_rates,u_temp,index] = ...
    HD_g_load_data_files(N,strain_spec,CK,0);
%% Fitting in steps: strate*exp(Q/RT) = A(sinh(alpha.sigma))^n
%# Loop over temp for estimating n
%# For that, the stress curves should be sorted for deformation temperature.
[temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx'],1);
%# The loop over temp
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
        t_X = log(sinh(ave_alpha(i)*t_stress(:,i)));
        t_Y = log(t_rate);
        curves_n(:,j) = [temp(l_temp,1) polyfit(t_X,t_Y,1) ...
            t_X' t_Y'];
    end
    
    %# Measure all data
    for i=1:N
        t_X = log(sinh(ave_alpha(i)*t_stress(:,i)));
        t_Y = log(t_rate);
        p_n(:,i) = polyfit(t_X,t_Y,1);
    end
    n(j,:) = p_n(1,:);
    
    l_temp = l_temp + num;
    j = j+1;
    if cond_temp(end)==1
        clear p_n
        l_temp = 0;
    end
    clear p_b
end
temp = temp(:,1);
if averaging_n
    ave_n = ones(1,size(n,2)) * mean(n(:));
else
    ave_n = mean(n);
end
clear t_temp l_temp t_rate t_stress num i j cond_temp t_cond_temp t_X p_var
%# Loop over strate for estimating Q
%# The stress curves should be sorted for deformation strain rate.
[temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx],2);

%# The loop over strate
[~,~,strate(:,2)] = unique(strate);
l_strate = 1;
j = 1;
Q_normaler = 1000;
while l_strate
    cond_strate = strate(:,2) == strate(l_strate,2);
    t_cond_strate = cond_strate;
    t_cond_strate(t_cond_strate==0) = [];
    num = size(t_cond_strate,1);

    t_temp = temp(l_strate:l_strate+num-1,1)';
    t_stress = ss_approx(l_strate:l_strate+num-1,:);
    
    %# Collect data for figure at the selected strain
    for i=index
        t_X = Q_normaler./(R*t_temp);
        t_Y = ave_n(i).*log(sinh(ave_alpha(i)*t_stress(:,i)));
        curves_Q(:,j) = [strate(l_strate,1) polyfit(t_X,t_Y',1) ...
            t_X t_Y'];
    end
    
    %# Measure all data
    for i=1:N
        t_X = Q_normaler./(R*t_temp);
        t_Y = ave_n(i).*log(sinh(ave_alpha(i)*t_stress(:,i)));
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

ave_n = ave_n';
%% Fit polynomial to averaged n,Q,lnA
p_ave_alpha = fit(strain,ave_alpha,p_order);
p_ave_n     = fit(strain,ave_n,p_order);
p_ave_Q     = fit(strain,ave_Q,p_order);
p_ave_lnA   = fit(strain,ave_lnA,p_order);

c_alpha= flip(coeffvalues(p_ave_alpha));
c_n    = flip(coeffvalues(p_ave_n));
c_Q    = flip(coeffvalues(p_ave_Q));
c_lnA  = flip(coeffvalues(p_ave_lnA)); 

poly_alpha  = feval(p_ave_alpha,strain);
poly_n      = feval(p_ave_n,strain);
poly_Q      = feval(p_ave_Q,strain);
poly_lnA    = feval(p_ave_lnA,strain);
%% Model results and errors
[mesh_strate,mesh_temp,mesh_n,mesh_Q,mesh_lnA,Z] = ...
    HD_g_mesher(strain,strate,temp,poly_n,poly_Q,poly_lnA,R,1000);
mesh_alpha = repmat(poly_alpha',size(temp));
sig = 1./mesh_alpha.*asinh((Z./exp(mesh_lnA)).^(1./mesh_n));
Correlation = sortrows([ss_approx(:),sig(:)]);
stress_m = Correlation(:,1);
stress_p = Correlation(:,2);
errors = HD_g_error_cal(stress_m,stress_p);
%% save workspace
if saving_f
    save('MAT_HD_sinh_conventional_fit',...
        'ss_approx','strain','strate','temp','u_rates','u_temp','index',...
        'strain_spec','averaging','p_order',...
        'ave_alpha','ave_n','ave_Q','ave_lnA','n','Q','lnA',...
        'poly_alpha','poly_n','poly_Q','poly_lnA',...
        'c_alpha','c_n','c_Q','c_lnA','curves_n','curves_Q',...
        'mesh_strate','mesh_temp','mesh_Q','mesh_lnA','mesh_alpha',...
        'stress_m','stress_p','sig','errors');
end
end