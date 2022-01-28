%-------------------------------------------------------------------------%
% Filename: figs_3_4_5_6_plot.m 
% Author: Juan M. Cardenas and Ben Adcock
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description: generates the plots for Figure 3, 4, 5, and 6.
%
% Inputs: 
% median_opt - option to plot median(1) or mean(0)
% semilog_opt - option to plot semilogy(1) or loglog(0)
% f - type of function (from 1 to 4) 
%
% Outputs: 
% figures 3,4,5,6
%-------------------------------------------------------------------------%

clear all; close all; clc;
addpath(genpath('../utils'))

for fig_num = 3:6

    median_opt   = 0;          % 1: Median   | 0: Mean
    semilog_opt  = 0;          % 1: semilogy | 0: loglog
    methods_fig  = [4 3 2];    % number of methods per figure

    if fig_num == 3
        f = 1;               
        init_row  = 1;       final_row = 5;

    elseif fig_num == 4
        f = 2;               
        init_row  = 1;       final_row = 4;

    elseif fig_num == 5
        f = 3;
        init_row  = 1;       final_row = 4;       

    elseif fig_num == 6
        f = 4;
        init_row  = 1;       final_row = 5;

    end
    
    file_name = ['../data/fig' ' ' sprintf('%d',fig_num) '/'  sprintf('Data_GDAS_fig_%d.mat',fig_num)];
    load(file_name) 
   
    % Recover data plot
    recover_GDAS_plot_data
   
    for row_num = init_row:final_row
        for col_num = 1:3
            
            if col_num == 1
                x_values    = Feval_matrix;
                y_values    = Error_matrix;
                num_methods = methods_fig(1);
            
            elseif col_num == 2  
                x_values    = Feval_vol_matrix;
                y_values    = Vol_matrix;
                num_methods = methods_fig(2);
                
            elseif col_num == 3
                x_values    = M_values_matrix;
                y_values    = Rate_matrix;
                num_methods = methods_fig(3);
            end
            
            fig = figure();
            
            for k = 1:num_methods
                
                plot(x_values(:,k,row_num),y_values(:,k,row_num),markers{13+k},...
                            'markersize',ms,...
                            'MarkerFaceColor',default_color(k,:),...
                            'MarkerEdgeColor',default_color(k,:),...
                            'LineWidth',lw);       
                hold on 
            end
            
            hold off 

            set_fonts_GDAS
            
            namefig    = sprintf('GDAS_fig_%d_%d_%d',fig_num,row_num,col_num);        
            foldername = ['../figs/fig' ' ' sprintf('%d',fig_num) '/'];
            saveas(fig, fullfile(foldername, namefig),'epsc');           
        end
    end
    
   % clear data
   clear all; close all; clc;
end