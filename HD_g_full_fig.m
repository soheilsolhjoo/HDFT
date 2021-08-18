% HD_g_full_fig(figure_reqest,tempC,fitting_figs,stress-srain_figs,corrlation_figs)
% 
% with a defaul call of: HD_g_full_fig([1,1,1,1,1],1,1,1,1);
%
% HD_full_fig is a general function of the Hot Deformation (HD) code. This
% function generates various figures related to the previously calibrated
% models.
%
% figure_reqest is in the format of [1x5] array of 0 or 1, with its
% elements respectively corresponding to:
% 1: power law model
% 2: exponential model
% 3: sine hyperbolic ST model (Sellars and Tegart)
% 4: revised ST model
% 5: trained ANN model
%
% Requesting figures for any of these sets work if the relevant files (with
% the auto saved names) are accessible.
% * tempC:              whether temperatures are in Kelvin(0) or Celcius(1)
% * fitting_figs:       generates plots used for fitting purposes
% * stress-srain_figs:  generates stress-strain curves
% * corrlation_fig:     generates the figure of correlation between the 
%                       measured and predicted stresses
% 
% For the ANN model, only the "corrlation_fig" is available.
%
%% DISCLAIMER
% This program is provided as is for free use.
%
%           Soheil Solhjoo
%           February 17, 2021
% -------------------------------------------------------------------------
function HD_g_full_fig(varargin)
idx = ~cellfun('isempty',varargin);
Defaults = {[1,1,1,1,1],1,1,1,1};
Defaults(idx) = varargin(idx);
[fig_req,tempC,fit_figs,ss_figs,corr_fig] = Defaults{:};

sym     = ['o','s','d','^','v','>','<','+','*','x','p','h','.'];
sym_corr= ["b^","ms","r+","go","cd"];

files = [...
    "MAT_HD_power_law_fit.mat",...
    "MAT_HD_exponential_fit.mat",...
    "MAT_HD_sinh_conventional_fit.mat",...
    "MAT_HD_sinh_revisited_fit.mat",...
    "MAT_HD_ANN.mat"];

for ml = 1:5 %main_loop
    if fig_req(ml) && isfile(files(ml))
        load(files(ml))
        max_S = max(stress_m)+min(stress_m);
        if ml ~= 5
            if averaging(1)
                mode = "const";
            else
                mode = "$f(\varepsilon)$";
            end
        else
            title_gen = join(["ANN model (",net.trainFcn, "), L=", string(L_num), ", N/L=", string(n_num)]);
        end
        if ml == 1
            
            title_gen = join(["Power-law model, $n'$: ",mode]);
            title_spec = join(["Power-law model, $\varepsilon = ",string(strain_spec),"$"]);
            
            data_t  = curves_n;
            data_r  = flip(curves_Q,2);
            xlabel_data_t = "$\ln(\sigma)$";
            ylabel_data_t = "$\ln(\dot{\varepsilon})$";
            ylabel_data_r = "$n'\ln(\sigma)$";
            
            dl_figs   = 3;
            dl_data   = ["n","Q","lnA"];
            dl_polys  = [poly_n,poly_Q,poly_lnA];
            dl_aves   = [ave_n,ave_Q,ave_lnA];
            dl_ylabs  = ["$n'$","$Q' (\rm kJ\cdot mol^{-1})$","$\ln(A^{'})$","$Q^{'}$"];
        elseif ml == 2
            title_gen = join(["Exponential model, $\beta$: ",mode]);
            title_spec = join(["Exponential model, $\varepsilon = ",string(strain_spec),"$"]);
            
            data_t  = curves_b;
            data_r  = flip(curves_Q,2);
            xlabel_data_t = "$\sigma$";
            ylabel_data_t = "$\ln(\dot{\varepsilon})$";
            ylabel_data_r = "$\beta\sigma$";
            
            dl_figs   = 3;
            dl_data   = ["b","Q","lnA"];
            dl_polys  = [poly_b,poly_Q,poly_lnA];
            dl_aves   = [ave_b,ave_Q,ave_lnA];
            dl_ylabs  = ["$\beta$","$Q'' (\rm kJ\cdot mol^{-1})$","$\ln(A^{''})$","$Q^{''}$"];
        elseif ml == 3
            if averaging(2)
                mode_n = "const";
            else
                mode_n = "$f(\varepsilon)$";
            end
            title_gen = join(["ST model, $\alpha$: ",mode,", $n:$ ",mode_n]);
            title_spec = join(["ST model, $\varepsilon = ",string(strain_spec),"$"]);
            
            data_t  = curves_n;
            data_r  = flip(curves_Q,2);
            xlabel_data_t = "$\ln(\sinh(\alpha\sigma))$";
            ylabel_data_t = "$\ln(\dot{\varepsilon})$";
            ylabel_data_r = "$n\ln(\sinh(\alpha\sigma)$";
            
            dl_figs   = 3;
            dl_data   = ["n","Q","lnA"];
            dl_polys  = [poly_n,poly_Q,poly_lnA];
            dl_aves   = [ave_n,ave_Q,ave_lnA];
            dl_ylabs  = ["$n$","$Q (\rm kJ\cdot mol^{-1})$","$\ln(A)$","$Q$"];
            
            if fit_figs
                figure; plot(strain,ave_alpha,'k-','LineWidth',1);
                xlabel("$\varepsilon$",'Interpreter','latex')
                ylabel("$\alpha$",'Interpreter','latex')
                title(title_gen,'Interpreter','latex');
                HD_g_figstyle
            end
        elseif ml == 4
            if averaging(2)
                mode_n = "const";
            else
                mode_n = "$f(\varepsilon)$";
            end
            title_gen = join(["R-ST model, method:", string(AQ),", $\alpha$: ",mode,", $n:$ ",mode_n]);
            title_spec = join(["R-ST model, method:", string(AQ),", $\varepsilon = ",string(strain_spec),"$"]);
            
            label_temp = label_temp_finder(u_temp',tempC);
        end
        if ~ismember(ml,[4,5]) && fit_figs
            %     ------------
            % linear plots for n'
            label_temp = label_temp_finder(data_t,tempC);
            HD_g_fig_generator(data_t,label_temp,'northwest',0);
            xlabel(xlabel_data_t,'Interpreter','latex')
            ylabel(ylabel_data_t,'Interpreter','latex')
            title(title_spec,'Interpreter','latex');
            HD_g_figstyle
            %     ------------
            % linear plots for Q'
            label = data_r(1,:);
            for i=1:size(data_r,2)
                label_rate(i) = {join([num2str(label(i)) 's^{-1}'])};
            end
            HD_g_fig_generator(data_r,label_rate,'northwest',0);
            xlabel_data_r = "$1000/RT$";
            xlabel(xlabel_data_r,'Interpreter','latex')
            ylabel(ylabel_data_r,'Interpreter','latex')
            title(title_spec,'Interpreter','latex');
            HD_g_figstyle
            
            dl_labs   = ["label_temp","label_rate","label_rate"];
            for k = 1:dl_figs
                dot_line(strain,eval(dl_data(k)),sym,eval(dl_labs(k)),dl_polys(:,k),dl_aves(:,k),title_gen,dl_ylabs(k))
            end
            ave_QlnA(strain,ave_Q,ave_lnA,poly_Q,poly_lnA,title_gen,dl_ylabs(4),dl_ylabs(3))
        end
        if ss_figs
            if ~exist('label_temp','var')
                label_temp = label_temp_finder(u_temp',tempC);
            end
            ss_fig(ss_approx,sig,mesh_strate,u_rates,u_temp,strain,sym,title_gen,label_temp)
        end
        if corr_fig
            corr(stress_m,stress_p,max_S,sym_corr(ml),title_gen,errors)
        end
    end
    clearvars -except tempC sym sym_corr fig_req files corr_fig fit_figs ss_figs
end
end
%% Complementary Functions
% stress-strain curves
function ss_fig(ss_approx,sig,mesh_strate,u_rates,u_temp,strain,sym,title_gen,label_temp)
ss_approx = sortrows([mesh_strate(:,1) ss_approx],1);
sig = sortrows([mesh_strate(:,1) sig],1);
ss_approx(:,1) = [];
sig(:,1) = [];
alpha = 0.75;
for i = 1:size(u_rates,1)
    set(0,'DefaultLegendAutoUpdate','off')
    figure; hold on
    for j=1:size(u_temp,1)
        l =  j + (i-1)*size(u_temp,1);
        scatter(strain,ss_approx(l,:),50,sym(j),...
            'DisplayName',string(label_temp(j)),'MarkerEdgeAlpha',alpha)
    end
    legend('Location','northeast');
    for j=1:size(u_temp,1)
        l =  j + (i-1)*size(u_temp,1);
        plot(strain,sig(l,:),'-k','LineWidth',1.5);
    end
    xlabel("$\varepsilon$",'Interpreter','latex')
    ylabel("$\sigma$(MPa)",'Interpreter','latex')
    title(join([title_gen,", $\dot{\varepsilon} = ",string(u_rates(i)),"$"]),'Interpreter','latex');
    hold off
    HD_g_figstyle
end
end

function dot_line(strain,n,sym,label_temp,poly_n,ave_n,fig_title,y_label)
figure; hold on
for i=1:size(n,1)
    scatter(strain,n(i,:),10,sym(i),...
        'DisplayName',string(label_temp(i)));
end
h = plot(strain,poly_n,'k-o','LineWidth',1,'DisplayName','Average');
legend('Location','northeast');%,'NumColumns',3);
delete(h)
plot(strain,poly_n,'k-','LineWidth',1);
plot(strain,ave_n,'ko','LineWidth',1);
xlabel("$\varepsilon$",'Interpreter','latex')
ylabel(y_label,'Interpreter','latex')
title(fig_title,'Interpreter','latex');
HD_g_figstyle
end

function corr(stress_m,stress_p,max_S,sym_corr,fig_title,errors)
figure; hold on
xlim([0 max_S])
ylim([0 max_S])
scatter(stress_m,stress_p,sym_corr);
xlabel("$\sigma_\textrm{m}\textrm{(MPa)}$",'Interpreter','latex')
ylabel("$\sigma_\textrm{p}\textrm{(MPa)}$",'Interpreter','latex')
title(fig_title,'Interpreter','latex');
X=[0 max_S]; Y=X;
plot(X,Y,'k-');
hold off
str_1 = join(["\textit{MARE}: ",num2str(round(errors.MARE,3))]);
str_2 = join(["$r^2$: ",num2str(round(errors.r2,3))]);
str = [str_1;str_2];
text(min(stress_m),0.8*max(stress_m),str,'Interpreter','latex','FontSize',16)
HD_g_figstyle
end

function ave_QlnA(strain,ave_Q,ave_lnA,poly_Q,poly_lnA,title_spec,QL,lnAL)
figure; hold on
title(title_spec,'Interpreter','latex');
xlabel("$\varepsilon$",'Interpreter','latex')
yyaxis left
ylabel(join([QL," $(k \rm J\cdot mol^{-1})$"]),'Interpreter','latex')
scatter(strain,ave_Q,20,'sb','DisplayName',QL);
yyaxis right
ylabel(lnAL,'Interpreter','latex')
scatter(strain,ave_lnA,20,'^r','DisplayName',lnAL);
legend('Interpreter','latex')
yyaxis left
plot(strain,poly_Q,'b')
yyaxis right
plot(strain,poly_lnA,'r')
HD_g_figstyle
end

function label_temp = label_temp_finder(data_t,tempC)
label = data_t(1,:);
for i=1:size(data_t,2)
    if tempC
        label_temp(i) = {[num2str(label(i)-273.15) ' °C']};
    else
        label_temp(i) = {[num2str(label(i)) ' K']};
    end
end
end