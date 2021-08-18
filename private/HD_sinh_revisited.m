% HD_sinh_revisited(CK,N,S,A,PO,AQ,norm,sf,mg):
% * CK:     Temperature Unit: 1 -> Celcius, 0 -> Kelvin
% * N:      Number of strain data points to be used for segmenting the
%           stress-strain curves
% * S:      Specified strain for collecting data for, e.g., plotting
% * A:      Averaging flag for parameters [alpha,n]: 1 -> yes, 0 -> no
% * PO:     Polynomial Order: acceptable values are 'poly1' to 'poly9'.
%           This one is for defining functions of strain.
% * AQ:     Following the suggested improving methods:
%           1 -> method 1 and  2 -> method 2
% * norm:   norm (Q_normalizer) is to normalize relevant data in
%           calculations of Q, which can be important in plotting data
% * sf:     saving flag. 1 -> yes, 0 -> no
% * mg:     model generator. 1 -> yes, 0 -> no
%
% For example: HD_sinh_revisited(1,100,0.05,[0,0],'poly9',1,1000,1,0)
% and the example show the default values for the function's arguments.
% 
% HD_sinh_revisited reads a hot deformation dataset and fit a revised
% version of the Sellars and Tegart's hyperbolic sine function on it. It 
% saves all calculations in 'MAT_HD_sinh_revisited.mat'.
%
% strate*exp(Q/RT) = A*(sinh(alpha*sigma))^n
% 
% The output includes calculated parameters and their approximations by
% polynomials. The stress values are also available. Comparisons can be
% made for:
% (measured)<-> (predicted)
% Q         <-> mesh_Q
% lnA       <-> mesh_lnA
% stress_m  <-> stress_p
% 
% Moreover, poly_Q and poly_lnA are strtuctures containing the coefficients
% of the fitted polynomials, with ".e", ".r" and ".t" corresponding to
% "strain", "strain-rate" and "temperature", respectively.
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
% -HD_sinh_conventional.m
% 
% Moreover, if A=[1,1] and mg=1, this function generates a self-standing
% function named 'M_HD_sinh_revisited_fit_model.m' based on the fitted
% model, which can be called as:
% 
% stress = ...
%       HD_sinh_revisited_fit_model(strain,strain_rate[1/s],temperature[K])
% 
%% DISCLAIMER
% This program is provided as is for free use.
%
%           Soheil Solhjoo
%           February 17, 2021
% -------------------------------------------------------------------------
function HD_sinh_revisited(varargin)
idx = ~cellfun('isempty',varargin);
Defaults = {1,100,0.05,[0,0],'poly9',1,1000,1,0};
Defaults(idx) = varargin(idx);
[CK,N,strain_spec,averaging,p_order,AQ,Q_normalizer,saving_f,model_gen] = ...
    Defaults{:};
R = 8.314;  %Gas constant

[averaging_alpha,averaging_n] = feval(@(x) x{:}, num2cell(averaging));
%% Load *.data files
[ss_approx,strain,strate,temp,u_rates,u_temp,index] = ...
    HD_g_load_data_files(N,strain_spec,CK,0);
%% Find | calculate ave_alpha and ave_n
[mesh_alpha,mesh_n,Q,lnA] = HD_sinh_conventional(CK,100,0.05,[0,0],'poly9',0);
if averaging_alpha
    mesh_alpha = mean(mesh_alpha(:));
end
if averaging_n
    mesh_n = mean(mesh_n(:));
end

%# Prepare data
[temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx'],1);
[temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx],2);
%% Fitting in steps: strate*exp(Q/RT) = Z = A (sinh(alpha * sigma))^n
%# This should be done in 6-steps
%# The first four steps are done by calling other functions
%# 5 & 6- find A and Q
%# NOTE: the order of steps 5 and 6 should be decided as an input:
%# AQ = 1: (5:Q & 6:A)
%# AQ = 2: (5:A & 6:Q)
%# 
%# Post-processing of the first estimations of A and Q
%# The first step of the post-processing is to define the steps 5 and 6, and
%# its general overview is as follows:
%# (a) post-processing of paramter in step 5
%# 1- the parameter in step 5 is averaged over strain rate
%# 2 & 3- the difference between the actual value and the average is
%# calculated, and approximated with a polynomial degree PO2 at each strain
%# 4- all constants are approximated with a polynomial degree "poly_degree"
%# for the entire range of strain
%# (b) finding the parameter in step 6
%# 1- estimate the parameter in step 6 purely as the proportionality
%# coefficient
%# 2- find its average over temperature for each set of strain rate (AQ=1)
%# or temperature (AQ=2);
%# 3- approximate the difference between the averages with a polynomial
%# degree 2 (1st dependence)
%# 4- find the difference between the average and the corresponding actual
%# values, and approximate it with a polynomial degree PO2 (2nd dependence)
if AQ == 1
    s5 = Q;
elseif AQ == 2
    s5 = lnA;
else
    error('Input:wrongdata', 'AQ must be either 1 or 2.');
end
%%%%%%%%
% step 5
s5_ave   = mean(s5);
s5_fit   = fit(strain,s5_ave',p_order);
s5_poly  = feval(s5_fit,strain);
s5_delta = s5 - s5_ave;
t_X = log(u_rates);
for i = 1:N
    t_Y = s5_delta(:,i);
    s5_delta_fit(:,i) = polyfit(t_X,t_Y,2);
end
s5_delta_fit_p0 = fit(strain,s5_delta_fit(3,:)',p_order);
s5_delta_fit_p1 = fit(strain,s5_delta_fit(2,:)',p_order);
s5_delta_fit_p2 = fit(strain,s5_delta_fit(1,:)',p_order);
t_X = log(strate);
s5_delta_poly = ...
    feval(s5_delta_fit_p2,strain)'.*t_X.^2+...
    feval(s5_delta_fit_p1,strain)'.*t_X+...
    feval(s5_delta_fit_p0,strain)';
mesh_s5 = s5_poly' + s5_delta_poly;
%%%%%%%%
% step 6
% In this step, the parameter is going to be treated purely as a
% proportionality constant; therefore, the meshgrids of strain, strin rate
% and temperature are required.
if AQ == 1
    [~,mesh_strate] = meshgrid(strain',strate(:,1)');
    [~,mesh_temp] = meshgrid(strain',temp');
    
    mesh_Q = mesh_s5;
    Z   = mesh_strate.*exp(Q_normalizer*mesh_Q/R./mesh_temp);
    D   = (sinh(mesh_alpha.*ss_approx)).^mesh_n;
    lnA = log(Z./D);
    
    % The loop over strain rate
    [~,~,strate(:,2)] = unique(strate);
    l_strate = 1;
    j = 1;
    lnA_delta = 0;
    while l_strate
        cond_strate = strate(:,2) == strate(l_strate,2);
        t_cond_strate = cond_strate;
        t_cond_strate(t_cond_strate==0) = [];
        num = size(t_cond_strate,1);
        
        lnA_tempo = lnA(cond_strate,:);
        lnA_ave_tempo(j,:) = mean(lnA_tempo);
        lnA_delta = lnA_delta + lnA_tempo - lnA_ave_tempo(j,:);
        
        l_strate = l_strate + num;
        j = j+1;
        if cond_strate(end)==1
            l_strate = 0;
        end
    end
    strate = strate(:,1);
    % Find difference between the averaged lnA and the whole set
    lnA_delta = lnA_delta / (j-1);
    [mesh_lnA,s6_ave_fit,...
        s6_ave_delta_fit_r_p0,s6_ave_delta_fit_r_p1,s6_ave_delta_fit_r_p2,...
        s6_ave_delta_fit_t_p0,s6_ave_delta_fit_t_p1,s6_ave_delta_fit_t_p2]=...
        s6_estimate(lnA_delta,lnA_ave_tempo,strain,N,ss_approx, ...
        p_order,u_rates,u_temp,strate,temp,AQ);
else
    [temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx],1);
    [~,mesh_strate] = meshgrid(strain',strate(:,1)');
    [~,mesh_temp] = meshgrid(strain',temp');
    
    mesh_lnA = mesh_s5;
    Z  = log(exp(mesh_lnA).*(sinh(mesh_alpha.*ss_approx)).^mesh_n./mesh_strate);
    C  = R.*mesh_temp/Q_normalizer;
    Q  = Z.*C;
    
    % The loop over temp
    [~,~,temp(:,2)] = unique(temp);
    l_temp = 1;
    j = 1;
    Q_delta = 0;
    while l_temp
        cond_temp = temp(:,2) == temp(l_temp,2);
        t_cond_temp = cond_temp;
        t_cond_temp(t_cond_temp==0) = [];
        num = size(t_cond_temp,1);
        
        Q_tempo = Q(cond_temp,:);
        Q_ave_tempo(j,:) = mean(Q_tempo);
        Q_delta = Q_delta + Q_tempo - Q_ave_tempo(j,:);
        
        l_temp = l_temp + num;
        j = j+1;
        if cond_temp(end)==1
            l_temp = 0;
        end
    end
    temp = temp(:,1);
    % Find difference between the averaged Q and the whole set
    Q_delta = Q_delta / (j-1);
    [mesh_Q,s6_ave_fit,...
        s6_ave_delta_fit_r_p0,s6_ave_delta_fit_r_p1,s6_ave_delta_fit_r_p2,...
        s6_ave_delta_fit_t_p0,s6_ave_delta_fit_t_p1,s6_ave_delta_fit_t_p2]=...
        s6_estimate(Q_delta,Q_ave_tempo,strain,N,ss_approx, ...
        p_order,u_rates,u_temp,strate,temp,AQ);
end
%% Save the fitted polynomials
if AQ == 1
    poly_Q.e        = s5_fit;
    poly_Q.r.p0     = s5_delta_fit_p0;
    poly_Q.r.p1     = s5_delta_fit_p1;
    poly_Q.r.p2     = s5_delta_fit_p2;
    poly_lnA.e      = s6_ave_fit;
    poly_lnA.r.p0   = s6_ave_delta_fit_r_p0;
    poly_lnA.r.p1   = s6_ave_delta_fit_r_p1;
    poly_lnA.r.p2   = s6_ave_delta_fit_r_p2;
    poly_lnA.t.p0   = s6_ave_delta_fit_t_p0;
    poly_lnA.t.p1   = s6_ave_delta_fit_t_p1;
    poly_lnA.t.p2   = s6_ave_delta_fit_t_p2;
else
    poly_lnA.e      = s5_fit;
    poly_lnA.r.p0   = s5_delta_fit_p0;
    poly_lnA.r.p1   = s5_delta_fit_p1;
    poly_lnA.r.p2   = s5_delta_fit_p2;
    poly_Q.e        = s6_ave_fit;
    poly_Q.r.p0     = s6_ave_delta_fit_r_p0;
    poly_Q.r.p1     = s6_ave_delta_fit_r_p1;
    poly_Q.r.p2     = s6_ave_delta_fit_r_p2;
    poly_Q.t.p0     = s6_ave_delta_fit_t_p0;
    poly_Q.t.p1     = s6_ave_delta_fit_t_p1;
    poly_Q.t.p2     = s6_ave_delta_fit_t_p2;
end
%% Model results and errors
Z   = mesh_strate.*exp(Q_normalizer*mesh_Q/R./mesh_temp);
sig = 1./mesh_alpha.*asinh((Z./exp(mesh_lnA)).^(1./mesh_n));
Correlation = sortrows([ss_approx(:),sig(:)]);
stress_m = Correlation(:,1);
stress_p = Correlation(:,2);
errors = HD_g_error_cal(stress_m,stress_p);
%% Export the fitted model as a function
if averaging_alpha && averaging_n && model_gen
syms x strain_s strate_s temp_s
s5_func_e = poly2sym(sym(coeffvalues(s5_fit)));
s5_func_r = ...
    poly2sym(sym(coeffvalues(s5_delta_fit_p2))).*log(strate_s(:)).^2 + ...
    poly2sym(sym(coeffvalues(s5_delta_fit_p1))).*log(strate_s(:))    + ...
    poly2sym(sym(coeffvalues(s5_delta_fit_p0)));
s5_func   = s5_func_e + s5_func_r;
s5_func   = subs(s5_func,x,strain_s);

s6_func_e = poly2sym(sym(coeffvalues(s6_ave_fit)));
s6_func_r = ...
    poly2sym(sym(coeffvalues(s6_ave_delta_fit_r_p2))).*log(strate_s(:)).^2 + ...
    poly2sym(sym(coeffvalues(s6_ave_delta_fit_r_p1))).*log(strate_s(:))    + ...
    poly2sym(sym(coeffvalues(s6_ave_delta_fit_r_p0)));
s6_func_t = ...
    poly2sym(sym(coeffvalues(s6_ave_delta_fit_t_p2))).*temp_s.^2 + ...
    poly2sym(sym(coeffvalues(s6_ave_delta_fit_t_p1))).*temp_s    + ...
    poly2sym(sym(coeffvalues(s6_ave_delta_fit_t_p0)));
s6_func   = s6_func_e + s6_func_r + s6_func_t;
s6_func   = subs(s6_func,x,strain_s);

if AQ == 1
    ZH          = strate_.*exp(Q_normalizer*s5_func/R./temp_s);
    sigma_fit   = 1./mesh_alpha.*asinh((ZH./exp(s6_func)).^(1./mesh_n));
else
    ZH          = strate_s.*exp(Q_normalizer*s6_func/R./temp_s);
    sigma_fit   = 1./mesh_alpha.*asinh((ZH./exp(s5_func)).^(1./mesh_n));
end
% Export as a m-code
matlabFunction(sigma_fit,'file','M_HD_sinh_revisited_fit_model.m')
end
%% save workspace
if saving_f
    save('MAT_HD_sinh_revisited_fit','AQ',...
        'ss_approx','strain','strate','temp','u_rates','u_temp','index',...
        'strain_spec','averaging','p_order','Q_normalizer','mesh_strate',...
        'mesh_alpha','mesh_n','Q','lnA','poly_Q','poly_lnA','mesh_Q','mesh_lnA',...
        'stress_m','stress_p','sig','errors');
end
end
%% Complementary Function
% Find the mesh_s6
function [mesh_s6,s6_ave_fit,...
    s6_ave_delta_fit_r_p0,s6_ave_delta_fit_r_p1,s6_ave_delta_fit_r_p2,...
    s6_ave_delta_fit_t_p0,s6_ave_delta_fit_t_p1,s6_ave_delta_fit_t_p2]= ...
    s6_estimate(s6_delta,s6_ave_tempo,strain,N,ss_approx, ...
    poly_degree,u_rates,u_temp,strate,temp,AQ)

s6_delta_ave       = mean(s6_delta)';
s6_delta_ave_fit   = fit(strain,s6_delta_ave,poly_degree);
s6_delta_ave_poly  = feval(s6_delta_ave_fit,strain);
s6_delta_ave_delta = s6_delta - s6_delta_ave_poly';
s6_ave       = mean(s6_ave_tempo)';
s6_ave_fit   = fit(strain,s6_ave,poly_degree);
s6_ave_poly  = feval(s6_ave_fit,strain);
s6_ave_delta = s6_ave_tempo - s6_ave_poly';

if AQ == 1
    % strain rate dependence
    t_X = log(u_rates);
    for i=1:N
        t_Y = s6_ave_delta(:,i);
        s6_ave_delta_fit(:,i) = polyfit(t_X,t_Y,2);
    end
    s6_ave_delta_fit_r_p0 = fit(strain,s6_ave_delta_fit(3,:)',poly_degree);
    s6_ave_delta_fit_r_p1 = fit(strain,s6_ave_delta_fit(2,:)',poly_degree);
    s6_ave_delta_fit_r_p2 = fit(strain,s6_ave_delta_fit(1,:)',poly_degree);
    t_X = log(strate);
    mesh_s6 = s6_ave_poly'.*ones(size(ss_approx));
    s6_ave_delta_poly = ...
        feval(s6_ave_delta_fit_r_p2,strain)'.*t_X.^2+...
        feval(s6_ave_delta_fit_r_p1,strain)'.*t_X+...
        feval(s6_ave_delta_fit_r_p0,strain)';
    mesh_s6 = mesh_s6 + s6_ave_delta_poly;
    % temperature dependence
    t_X = u_temp;
    for i=1:N
        t_Y = s6_delta_ave_delta(:,i);
        s6_ave_delta_fit(:,i) = polyfit(t_X,t_Y,2);
    end
    s6_ave_delta_fit_t_p0 = fit(strain,s6_ave_delta_fit(3,:)',poly_degree);
    s6_ave_delta_fit_t_p1 = fit(strain,s6_ave_delta_fit(2,:)',poly_degree);
    s6_ave_delta_fit_t_p2 = fit(strain,s6_ave_delta_fit(1,:)',poly_degree);
    t_X = temp;
    s6_ave_delta_poly = ...
        feval(s6_ave_delta_fit_t_p2,strain)'.*t_X.^2+...
        feval(s6_ave_delta_fit_t_p1,strain)'.*t_X+...
        feval(s6_ave_delta_fit_t_p0,strain)';
    mesh_s6 = mesh_s6 + s6_ave_delta_poly;
else
    % temperature dependence
    t_X = u_temp;
    for i=1:N
        t_Y = s6_ave_delta(:,i);
        s6_ave_delta_fit(:,i) = polyfit(t_X,t_Y,2);
    end
    s6_ave_delta_fit_t_p0 = fit(strain,s6_ave_delta_fit(3,:)',poly_degree);
    s6_ave_delta_fit_t_p1 = fit(strain,s6_ave_delta_fit(2,:)',poly_degree);
    s6_ave_delta_fit_t_p2 = fit(strain,s6_ave_delta_fit(1,:)',poly_degree);
    t_X = temp;
    mesh_s6 = s6_ave_poly'.*ones(size(ss_approx));
    s6_ave_delta_poly = ...
        feval(s6_ave_delta_fit_t_p2,strain)'.*t_X.^2+...
        feval(s6_ave_delta_fit_t_p1,strain)'.*t_X+...
        feval(s6_ave_delta_fit_t_p0,strain)';
    mesh_s6 = mesh_s6 + s6_ave_delta_poly;
    
    % strain rate dependence
    t_X = log(u_rates);
    for i=1:N
        t_Y = s6_delta_ave_delta(:,i);
        s6_ave_delta_fit(:,i) = polyfit(t_X,t_Y,2);
    end
    s6_ave_delta_fit_r_p0 = fit(strain,s6_ave_delta_fit(3,:)',poly_degree);
    s6_ave_delta_fit_r_p1 = fit(strain,s6_ave_delta_fit(2,:)',poly_degree);
    s6_ave_delta_fit_r_p2 = fit(strain,s6_ave_delta_fit(1,:)',poly_degree);
    t_X = log(strate);
    s6_ave_delta_poly = ...
        feval(s6_ave_delta_fit_r_p2,strain)'.*t_X.^2+...
        feval(s6_ave_delta_fit_r_p1,strain)'.*t_X+...
        feval(s6_ave_delta_fit_r_p0,strain)';
    mesh_s6 = mesh_s6 + s6_ave_delta_poly;
end
end