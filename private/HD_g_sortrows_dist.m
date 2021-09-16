% [temp,strate,stress] = HD_g_sortrows_dist([temperature,strain_rate,stress],row)
% 
% HD_g_sortrows_dist is a general function of HDFT. This function sorts 
% input data with three input matrices "inputMAT", with the format: 
% [temperature,strate,stress].
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
function [temp,strate,stress] = HD_g_sortrows_dist(inputMAT,row)
inputMAT= sortrows(inputMAT,row);%sort data for row
temp    = inputMAT(:,1);                %sorted temperature list
strate  = inputMAT(:,2);              %sorted strain rate list
stress  = inputMAT(:,3:end);          %sorted stress matrix
end