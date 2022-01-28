%-------------------------------------------------------------------------%
% Filename: figs_2_run.m 
% Author: Juan M. Cardenas and Ben Adcock
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description: generates the data and plots for Figure 2.
%
% Inputs: 
% median_opt - option to plot median(1) or mean(0)
% semilog_opt - option to plot semilogy(1) or loglog(0)
% f - type of function (from 1 to 4) 
%
% Outputs: 
% figure 2
%-------------------------------------------------------------------------% 
close all; clear all; clc;
addpath(genpath('../utils'))

% Generate grid 
fig_num    = 2;
K = 60000;
d = 2;
Z = 2*rand(K,d)-1;
indices_D = [1:K]';


for col_num = 1:3
    
    if col_num == 1   
        f = 2;   
        r_d = 7/10;
        domain = 1; tol = 0; a = 0; b = 0;
    
    elseif col_num == 2
        f = 1; 
        domain = 1; tol = 0; a = 0; b = 0;
        
    elseif col_num == 3
        f = 3; 
        domain = 2; tol = 0; a = 1.8/10; b = 0.72;
    end
    
    % function
    switch f
        case 1          
            % High-dim: Annulus
            Y      = sum(Z.^2,2);
            f_grid = log(8*Y)-(2)*Y;             
        case 2          
            % high-dim: Prism
            y1     = Z(:,1);    
            y2     = Z(:,2);
            f_grid = ((1/r_d)^2 - 1./(y1.^2 + y2.^2)).*exp(-(sum(Z,2)/(2*d)));
        case 3   
            ind    = [1:d] + 1;
            v_1    = (Z + (((-1).^ind)./ind)).^2;
            f_grid = prod( (d/4)./((d/4) + v_1),2);
        case 4          
            % High-dim: New f_1
            Y      = sum(Z.^2,2);
            f_grid = g_dim(d)*log((16/d)*Y)-(4/d)*Y;              
    end
    
    % Domain Omega
    if domain == 1
        I_Omega = indices_D(f_grid >= tol);
    else
        % domain 2 
        I_om1 = indices_D(a<= f_grid);
        I_om2 = indices_D(f_grid <= b);
        I_Omega = intersect(I_om1,I_om2);
    end
    K_Omega = length(I_Omega);
    
    fig = figure();
    x   = Z(I_Omega,1);
    y   = Z(I_Omega,2);
    shp = alphaShape(x,y,0.15);
    plot(shp,'EdgeColor','none')
    xlim([-1,1])
    ylim([-1,1])
 
	ax            = gca;   
	ax.FontSize   = 15;                 ax.LineWidth  = 3;
	ax.YGrid      = 'on';               ax.XGrid      = 'on';
	ax.XMinorTick = 'on';               ax.YMinorTick = 'off';
    ax.YTick = [-1 -0.5 0 0.5 1];
    ax.XTick = [-1 -0.5 0 0.5 1];
	
    hold off 
    
    beautify_plot
            
	namefig    = sprintf('GDAS_fig_%d_%d',fig_num,col_num);        
	foldername = ['../figs/fig' ' ' sprintf('%d',fig_num) '/'];
	saveas(fig, fullfile(foldername, namefig),'epsc');
    
end