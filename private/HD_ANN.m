% HD_ANN(CK,[L_num,n_num],algorithm,sf,mg):
% * CK:             Temperature Unit: 1 -> Celcius, 0 -> Kelvin
% * sf:             saving flag. 1 -> yes, 0 -> no
% * mg:             model generator. 1 -> yes, 0 -> no
% * [L_num,n_num]:  number of layers and neurons/layer
% * algorithm:      the algorithm for training the ANN. Further details can
%                   be found in MathWorks website:
%                   (1) https://mathworks.com/help/deeplearning/ug/train-and-apply-multilayer-neural-networks.html
%                   (2) https://mathworks.com/help/deeplearning/ref/feedforwardnet.html
%
% For example: HD_ANN(1,[2,5],'trainbr',1,1)
% and the example show the default values for the function's arguments.
%
% HD_ANN reads a hot deformation dataset and fit an ANN on it, and save all
% calculations in 'MAT_HD_ANN.mat'.
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
% -HD_g_error_cal.m
%
%% DISCLAIMER
%   HDFT  Copyright (C) 2021  Soheil Solhjoo
%   This program comes with ABSOLUTELY NO WARRANTY.
%   This is free software, and you are welcome to redistribute it under certain conditions.
%   Check "copyright.txt" in the main folder.
%
%           Soheil Solhjoo
%           February 17, 2021
% -------------------------------------------------------------------------
function HD_ANN(varargin)
%# Check the input argument
idx = ~cellfun('isempty',varargin);
Defaults = {1,[2,5],'trainbr',1,1};
Defaults(idx) = varargin(idx);
[CK,num,algorithm,saving_f,model_gen] = Defaults{:};
%# Distribute the assigned values
[L_num,n_num] = feval(@(x) x{:}, num2cell(num));
%% Input data
%# Loading in the raw data
data = HD_g_load_data_files(0,0,CK,1);
%# Formating the data for the NN
Temp = data.ss(:,1)';
deps = data.ss(:,2)';
eps  = data.ss(:,3)';
sig  = data.ss(:,4)';
N    = data.N;
%# Building data for NN
input  = [eps(:) deps(:) Temp(:)]';
stress_m = sig(:)';
%% Train the NN
net = fitnet(ones(1,L_num)*n_num,algorithm);
[net,tr] = train(net,input,stress_m);
%% Find the trained NN's approximations
stress_p = net(input);
%% Error Measurements
errors = HD_g_error_cal(stress_m,stress_p);
%% save the calculations
if saving_f
    
    % Collect data and make estimations for plotting stress-strain curves.
    [ss_approx,strain,strate,temp,u_rates,u_temp,~] = ...
        HD_g_load_data_files(floor(size(sig,2)/N),min(eps),CK,0);
    [temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx'],1);
    [temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx],2);
    [mesh_strain,mesh_strate] = meshgrid(strain',strate(:,1)');
    [~,mesh_temp]   = meshgrid(strain',temp');
    input_est = [mesh_strain(:) mesh_strate(:) mesh_temp(:)]';
    sig = reshape(net(input_est),size(ss_approx));
    
    % SAVE data
    save('MAT_HD_ANN','ss_approx','strain','u_rates','u_temp','mesh_strate','sig',...
        'net','tr','input','stress_m','stress_p','errors','L_num','n_num');
end
%% Export the fitted model as a function
if model_gen
    genFunction(net,'M_HD_ANN_fit_model.m');
end
end
