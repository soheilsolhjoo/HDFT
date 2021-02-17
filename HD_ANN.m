%% DISCLAIMER
% This program is provided as is for free use. 
% 
%           Soheil Solhjoo
%           FSE, University of Groningen, The Netherlands
%           February 17, 2021
%% Description
% HD_ANN reads a hot deformation dataset and fit an ANN on it, and save all
% calculations in 'trained_ANN.mat'.
% 
% The functions should be called as:
% HD_ANN(filename,[t,d,e,s],[L_num,n_num],algorithm):
% * CSV_filename:   the name of the file in a columnar format with each
%                   column assigned to either temperature (K),
%                   strain-rate (1/s), strain or stress (MPa).
% * [t,d,e,s]:      the corresponding collumns for temperature,
%                   strain-rate, strain and stress, respectively.
% * [L_num,n_num]:  number of layers and neurons/layer
% * algorithm:      the algorithm for training the ANN. Further details can
%                   be found in MathWorks website:
%                   (1) https://mathworks.com/help/deeplearning/ug/train-and-apply-multilayer-neural-networks.html
%                   (2) https://mathworks.com/help/deeplearning/ref/feedforwardnet.html
% 
% For example: HD_ANN('raw_data.csv',[1,2,3,4],[2,5],'trainbr')
% 
% NOTE: the two first arguments are mandatory. The default values for the
% last two ones are [L_num = 2, n_num = 5] and algorithm = 'trainbr'.
% -------------------------------------------------------------------------
function HD_ANN(varargin)%CSV_filename,tdes,num,algorithm)
%# Check the input argument
idx = ~cellfun('isempty',varargin);
if size(idx,2) < 2 || size(idx,2) > 4 || size(varargin{2},2) ~= 4
    error('Input:incorrectStruc',...
            'The input is invalid.');
end
Defaults = {[],[],[2,5],'trainbr'};
Defaults(idx) = varargin(idx);
[filename,tdes,num,algorithm] = Defaults{:};
%# Distribute the assigned values
[t,d,e,s] = feval(@(x) x{:}, num2cell(tdes));
[L_num,n_num] = feval(@(x) x{:}, num2cell(num));
%% Input data
%# Loading in the raw data
data = readmatrix(filename);
%# Formating the data for the NN
Temp = data(:,t);
deps = data(:,d);
eps  = data(:,e);
sig  = data(:,s);
%# Transposing the data
Temp = Temp.';
deps = deps.';
eps  = eps.';
sig  = sig.';
%# Building data for NN
input  = [eps(:) deps(:) Temp(:)]';
stress_m = sig(:)';
%% Train the NN
net = fitnet(ones(1,L_num)*n_num,algorithm);
[net,tr] = train(net,x,stress_m);
%% Find the trained NN's approximations
stress_p = net(input);
%% Error Measurements
% Mean Average Relative Error
MARE = 1/size(stress_m,2)*...
    sum(abs((stress_m-stress_p)./stress_m));
% Coefficient of determination
r2   = 1 - ...
    (sum((stress_m-stress_p).^2)) / ...
    (sum((stress_m-mean(stress_m)).^2));
% Pearson correlation coefficient
rxy  = sum((stress_p-mean(stress_p)).*(stress_m-mean(stress_m))) / ...
    (sqrt(sum((stress_p-mean(stress_p)).^2)) * ...
    sqrt(sum((stress_m-mean(stress_m)).^2)));
errors.MARE = MARE;
errors.r2 = r2;
errors.rxy = rxy;
%% save the calculations
save('trained_ANN','net','tr','input','stress_m','stress_p','errors');