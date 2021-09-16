% load_data_files is a general function of HDFT. This function reads 
% "*.data" files to be used in other HDFT files.
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
function [ss_approx,strain,strate,temp,u_rates,u_temp,index] = ...
    HD_g_load_data_files(N,strain_spec,CK,ANN)
% cent is to check whether the temperature is in Kelvin (0) or Celcius(1).
if CK
    temp_C = 273.15;
else
    temp_C = 0;
end
%# Prepare List of *.DATA files
files = dir('*.data');
names = {files.name};
names = split(names','.');
if isempty(names)
    error('Input:missingFiles', 'The *.data files are not found.');
end
names = str2double(names(:,1));
names = sort(names);
names(isnan(names(:)))=[];

%# Read stress-strain curves
ns = size(names,1);
temp = zeros(ns,1);
strate = zeros(ns,1);

%# Check if the function is called from within HD_ANN
if ANN
stress_strain = zeros(1,4);
    for i=1:ns
        input = dlmread([int2str(names(i)),'.data'], '\t');
        si = size(input,1) - 2;
        ss = size(stress_strain,1) + 1;
        stress_strain(ss:ss+si-1,3:4) = input(3:end,:);
        stress_strain(ss:ss+si-1,1) = input(1,1) + temp_C;
        stress_strain(ss:ss+si-1,2) = input(2,1);
    end
    ss_approx.ss = stress_strain(2:end,:);
    ss_approx.N  = ns;
    return
end

for i=1:ns
    input = dlmread([int2str(names(i)),'.data'], '\t');
    temp(i,1) = input(1,1) + temp_C;
    strate(i,1) = input(2,1);
    ss = sortrows(input(3:end,:));
%     stress_strain(:,:,i) = input(3:end,:);
    stress_strain{i} = ss;
    min_strain(i) = ss(1,1);
    max_strain(i) = ss(end,1);
end

%# Find unique values strate and temp
u_rates = unique(strate);
u_temp = unique(temp);

%# Find the smallest recorded strain and cut all curves for that
max_strain = min(max_strain);
min_strain = max(min_strain);
strain = linspace(0,max_strain,N)';
index = find(round(strain,3)==strain_spec);
ss_approx = zeros(N,ns);
for i=1:ns
    ss = stress_strain{i};
    ss(ss(:,1) < min_strain,:) = [];
    ss(ss(:,1) > max_strain,:) = [];
    ss_approx(:,i) = ...
        interp1(ss(:,1),ss(:,2),strain,'spline','extrap');
end
end