%-------------------------------------------------------------------------%
% Filename: figs_7_plot.m 
% Author: Juan M. Cardenas and Ben Adcock
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description: generates the plots for Figure 7.
%
% Inputs: 
% median_opt - option to plot median(1) or mean(0)
% semilog_opt - option to plot semilogy(1) or loglog(0)
% f - type of function (from 1 to 4) 
%
% Outputs: 
% figures 7
%-------------------------------------------------------------------------% 
clear all; close all; clc;
addpath(genpath('../utils'))

for col_num = 1:3
    
    fig_num =     7;
    
    median_opt   = 0;          % 1: Median   | 0: Mean
    semilog_opt  = 0;          % 1: semilogy | 0: loglog
    methods_fig  = [4 3 2];    % number of methods per figure
    
    if col_num == 1
       f = 1;      
    elseif col_num == 2 
       f = 2;    
    elseif col_num == 3                   
       f = 4;   
    end
    
    file_name = ['../data/fig' ' ' sprintf('%d',fig_num) '/' sprintf('Data_GDAS_fig_%d_f_%d_dim2.mat',fig_num,f)];
    load(file_name) 
   
    % Recover data plot
    recover_GDAS_plot_data
   
    for row_num = 1:5
             
        if row_num == 1
            x_values    = Feval_matrix;
            y_values    = Error_matrix;
            num_methods = methods_fig(1);
            
        elseif row_num == 2  
            x_values    = Feval_vol_matrix;
            y_values    = Vol_matrix;
            num_methods = methods_fig(2);
                
        elseif row_num == 3
            x_values    = M_values_matrix;
            y_values    = Rate_matrix;
            num_methods = methods_fig(3);
        
        elseif row_num == 4    
            x_values    = l_steps;
            y_values    = J_vals_matrix;
            num_methods = methods_fig(1);
            
        elseif row_num == 5
            x_values    = l_steps;
            y_values    = C_vals_matrix;
            num_methods = methods_fig(1);            
        end
            
        fig = figure();
            
        for k = 1:num_methods
     
            plot(x_values(:,k),y_values(:,k),markers{13+k},...
                 'markersize',ms,...
                 'MarkerFaceColor',default_color(k,:),...
                 'MarkerEdgeColor',default_color(k,:),...
                 'LineWidth',lw);       
            hold on 
            
        end
            
        hold off 

        set_fonts_GDAS_fig_7

        namefig    = sprintf('GDAS_fig_%d_%d_%d',fig_num,row_num,col_num);        
        foldername = ['../figs/fig' ' ' sprintf('%d',fig_num) '/'];          
             
        saveas(fig, fullfile(foldername, namefig),'epsc');           
    end
    
   % clear data
   clear all; close all; clc;
end