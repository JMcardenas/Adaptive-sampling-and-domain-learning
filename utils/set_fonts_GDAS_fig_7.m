%-------------------------------------------------------------------------%
% Filename: set_fonts_GDAS_fig_7.m
% Authors: Juan M. Cardenas  
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description:  This script set fonts fig_7
%
%-------------------------------------------------------------------------% 
 
if row_num == 1
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    h = legend('MC-LS','ASUD-LS','ASUD-ALS','ASGD-LS');
    
elseif row_num == 2
    set(gca,'xscale','log'); 
    h = legend('MC-LS','ASUD-LS','ASUD-ALS');
    
elseif row_num == 3
    h = legend('MC-LS','ASUD-LS');   
    
elseif row_num == 4 
    set(gca,'yscale','log');
    h = legend('MC-LS','ASUD-LS','ASUD-ALS','ASGD-LS');
    
elseif row_num == 5
    set(gca,'yscale','log');
    h = legend('MC-LS','ASUD-LS','ASUD-ALS','ASGD-LS');
    
end

set(h,'Interpreter','latex','location','best');  

ax            = gca;   
ax.FontSize   = 13;                 ax.LineWidth  = 1.5;
ax.YGrid      = 'on';               ax.XGrid      = 'on';
ax.XMinorTick = 'on';               ax.YMinorTick = 'off';

if semilog_opt == 1
    ax.XMinorGrid = 'off';       
    ax.YMinorGrid = 'off';
else 
    ax.XMinorGrid = 'on'; 
    ax.YMinorGrid = 'off';
end
