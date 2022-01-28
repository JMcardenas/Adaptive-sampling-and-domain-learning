%-------------------------------------------------------------------------%
% Filename: figs_3_4_5_6_run.m 
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description: generates all data for Figure 3,4,5,6.
%
% Inputs: 
% median_opt - option to plot median(1) or mean(0)
% semilog_opt - option to plot semilogy(1) or loglog(0)
% f - type of function (from 1 to 4) 
% save_file - option to save the file
% plt_domain - option to plot the domain
% N_max - maximum size of index set. 
% max_div - divisions that define the number of iterations
% Trials - number of trials 
% dim - dimension

% Outputs: 
% data for figure 3,4,5,6.
%-------------------------------------------------------------------------% 

clear all; close all; clc;

for fig_num = 3:6
  
   addpath(genpath('../utils'))
   
   Trials     = 50;
   index_type = 1;
   save_file  = 1;
   plt_domain = 1;
   N_max      = 2000;
   max_div    = 18;
   
   if fig_num == 3
       domain = 1;      tol = 0;             a   = 0;    b = 0;
       f      = 1;      dim = [2 3 5 10 15]; r_f = 2;
       
   elseif fig_num == 4
       domain = 1;      tol = 0;            a   = 0;     b = 0;
       f      = 2;      dim = [2 3 4 5];    r_f = 1;
       
   elseif fig_num == 5
       domain = 1;       tol = 0;           a   = 0;     b = 0;
       f      = 3;       dim = [2 3 4 5];   r_f = 1;
       
   elseif fig_num == 6
       domain = 2;       tol = 0;             a   = 1.8/10;      b = 0.72;                      
       f      = 4;       dim = [2 3 5 10 15]; r_f = 1;
   end
   
   % Run GDAS
   GDAS_data_figs
   
   % clear data
   clear all; close all; clc;
end 
 
 
 