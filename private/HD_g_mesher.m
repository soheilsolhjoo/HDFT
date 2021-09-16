% HD_g_mesher is a general function of HDFT. This function takes various 
% parameters as the input argument and create a representative mesh 
% (2D matrix) of them. It also calculates the Zenner-Holomon parameter.
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
function [mesh_strate,mesh_temp,mesh_c,mesh_Q,mesh_lnA,Z] = ...
    HD_g_mesher(strain,strate,temp,poly_c,poly_Q,poly_lnA,R,Q_normalizer)
[~,mesh_strate] = meshgrid(strain',strate(:,1)');
[~,mesh_temp]   = meshgrid(strain',temp');
mesh_c          = repmat(poly_c',size(temp));
mesh_Q          = repmat(poly_Q',size(temp));
mesh_lnA        = repmat(poly_lnA',size(temp));

Z   = mesh_strate.*exp(Q_normalizer*mesh_Q/R./mesh_temp);
end