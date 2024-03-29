% errors = HD_g_error_cal(Measured_Data,Predicted_Data)
% 
% error_cal is a general function of HDFT. This function calculates three
% different errors between measured data (x_m) and predicted values (x_p).
% The output is a structure with three components of Mean Average Relative 
% Error (MARE), coefficient of determination (r2), and Pearson correlation 
% coefficient (rxy).
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
function errors = HD_g_error_cal(x_m,x_p)
% Mean Average Relative Error
MARE = 1/length(x_m)*...
    sum(abs((x_m-x_p)./x_m));
% Coefficient of determination
r2   = 1 - ...
    (sum((x_m-x_p).^2)) / ...
    (sum((x_m-mean(x_m)).^2));
% Pearson correlation coefficient
rxy  = sum((x_p-mean(x_p)).*(x_m-mean(x_m))) / ...
    (sqrt(sum((x_p-mean(x_p)).^2)) * ...
    sqrt(sum((x_m-mean(x_m)).^2)));

errors.MARE = MARE;
errors.r2 = r2;
errors.rxy = rxy;
end