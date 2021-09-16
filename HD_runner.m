clear; close all;

% parameters
CK = 1;
N = 100;
picked_strain = 0.05;
const = 1;
poly_type = 'poly9';
revised_ST_method = 1;
Q_normliazer = 1000;
save_flag = 1;
model_gen = 1;
ANN_struct = [2,5];
ANN_alg = 'trainbr';

% functions call
HD_f_power_law(CK,N,picked_strain,const,poly_type,save_flag);
HD_f_exponential(CK,N,picked_strain,const,poly_type,save_flag);
HD_f_sinh_conventional(CK,N,picked_strain,[const,const],poly_type,save_flag);
HD_f_sinh_revisited(CK,N,picked_strain,[const,const],poly_type,revised_ST_method,Q_normliazer,save_flag,model_gen);
HD_f_ANN(CK,N,ANN_struct,ANN_alg,save_flag,model_gen);

% produce plots
HD_g_full_fig
