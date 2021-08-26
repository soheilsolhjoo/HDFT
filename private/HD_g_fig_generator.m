% HD_g_fig_generator is a general function of the Hot Deformation (HD)
% code. This function generates requested figures with (up to 13) various 
% data sets.
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
function f = HD_g_fig_generator(data,label,legend_position,latex_legend)
set(0,'DefaultLegendAutoUpdate','off')
sym     = ['o','s','d','^','v','>','<','+','*','x','p','h','.'];
polyc   = data(2:3,:);
data    = data(4:end,:);
list_s  = size(data,1)/2;
x       = data(1:list_s,:);
y       = data(list_s+1:end,:);
f = figure; hold on
for i=1:size(data,2)
    x_temp = x(:,i);
    y_temp = y(:,i);
    scatter(x_temp,y_temp,50,sym(i),'filled',...
        'DisplayName',string(label(i)))
end
if latex_legend
    legend('Location',legend_position,'Interpreter','latex');
else
    legend('Location',legend_position);
end
for i=1:size(data,2)
    x_temp = x(:,i);
    p_temp = polyc(:,i);
    plot(x_temp,polyval(p_temp,x_temp),'--k')
end
end
