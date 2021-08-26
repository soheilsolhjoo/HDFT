% HD_g_figstyle is a general function of the Hot Deformation (HD) code. 
% This function applies a predefined style on the generated figures.
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
function HD_g_figstyle
set(gca,'TickLabelInterpreter', 'latex','FontSize', 16);
ax = gca;
ax.TitleFontSizeMultiplier = 0.9;
axis xy

x0 = 100;
y0 = x0;
x1 = x0 +400;
y1 = y0 +300;
set(gcf, 'Position', [x0 y0 x1 y1])
set(gcf, 'PaperPositionMode', 'auto');

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2) + 0.01;
ax_width = outerpos(3) - ti(1) - ti(3) - 0.005;
ax_height = outerpos(4) - ti(2) - ti(4) -0.02;
ax.Position = [left bottom ax_width ax_height];
end
