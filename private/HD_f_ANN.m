% HD_f_ANN(CK,N,[L_num,n_num],algorithm,sf,mg):
% * CK:             Temperature Unit: 1 -> Celcius, 0 -> Kelvin
% * N:              Number of strain data points to be used for segmenting the stress-strain curves
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
% calculations in 'HDFT_fit.mat' as a stucture named "ANN#_#", with # being
% L_num and n_num, respectively.
% 
% If mg = 1, the trained network will be saved as function named
% "HD_ANN_fit_model.m".
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
%           v1: February 17, 2021
%           v2: September 16, 2021
% -------------------------------------------------------------------------
function HD_f_ANN(varargin)
%# Check the input argument
idx = ~cellfun('isempty',varargin);
Defaults = {1,100,[2,5],'trainbr',1,1};
Defaults(idx) = varargin(idx);
[CK,N,num,algorithm,saving_f,model_gen] = Defaults{:};
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
% N    = data.N;
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
        HD_g_load_data_files(N,min(eps),CK,0);
    [temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx'],1);
    [temp,strate,ss_approx] = HD_g_sortrows_dist([temp,strate,ss_approx],2);
    [mesh_strain,mesh_strate] = meshgrid(strain',strate(:,1)');
    [~,mesh_temp]   = meshgrid(strain',temp');
    input_est = [mesh_strain(:) mesh_strate(:) mesh_temp(:)]';
    sig = reshape(net(input_est),size(ss_approx));
    
    HDFT_fit_temp.type = 5;
    data = ["input,call,ANN,mesh,corr";
        "ss_approx,strain,strate,temp,u_rates,u_temp";
        "CK,N,L_num,n_num";
        "net,tr";
        "mesh_strate,mesh_temp";
        "stress_m,stress_p,sig";];
    var_name = strsplit(data(1),',');
    for i = 2:size(var_name,2)+1
        list = strsplit(data(i),',');
        for j = 1:size(list,2)
            HDFT_fit_temp.(var_name(i-1)).(list(j)) = eval(list(j));
        end
    end
    HDFT_fit_temp.errors = errors;
    input = HDFT_fit_temp.input;
    HDFT_fit_temp = rmfield(HDFT_fit_temp,'input');
    toSave = ["ANN" + num2str(L_num) + "_" + num2str(n_num)];
    HD_g_assign(toSave,HDFT_fit_temp)
    
    filename = "HDFT_fit.mat";
    if ~isfile(filename)
        save('HDFT_fit','input')
    end
    save('HDFT_fit',toSave,'-append')
end
%% Export the fitted model as a function
if model_gen
    genFunction(net,'HD_ANN_fit_model.m');
end
end