% [temp,strate,stress] = HD_g_sortrows_dist([temperature,strain_rate,stress],row)
% 
% HD_g_sortrows_dist is a general function of the Hot Deformation (HD)
% code. This function sorts input data with three input matrices
% "inputMAT", with the format: [temperature,strate,stress].
% 
%% DISCLAIMER
% This program is provided as is for free use.
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