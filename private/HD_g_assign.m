% HD_g_assign(var_name,value,extract)
% 
% assign is a general function of HDFT. This function generates a variable
% with name "var_name" in workspace of the caller function and assign
% "value" to it.
% 
% extract = 0 : to assign a value to a variable         HD_g_assign(var_name,value)
% extract = 1 : to extract all fields of a structure    HD_g_assign(~,structure,1)
% 
%% DISCLAIMER
%   HDFT  Copyright (C) 2021  Soheil Solhjoo
%   This program comes with ABSOLUTELY NO WARRANTY.
%   This is free software, and you are welcome to redistribute it under certain conditions.
%   Check "copyright.txt" in the main folder.
% 
%           Soheil Solhjoo
%           September 16, 2021
% -------------------------------------------------------------------------
function HD_g_assign(varargin)
%# Check the input argument
idx = ~cellfun('isempty',varargin);
Defaults = {[],[],0};
Defaults(idx) = varargin(idx);
[var_name,value,extract] = Defaults{:};

if ~extract
    assignin('caller', var_name, value);
else
    struct_fn = string(fieldnames(value));
    struct_fn(struct_fn(:)=="type") = [];
    for i=1:numel(struct_fn)
        data_temp = value.(struct_fn(i));
        if isstruct(data_temp)
            data_temp_fn = string(fieldnames(data_temp));
            for j=1:numel(data_temp_fn)
                assignin('caller',data_temp_fn(j),data_temp.(data_temp_fn(j)));
            end
        else
            assignin('caller',struct_fn(i),data_temp);
        end
    end
end
end