%-------------------------------------------------------------------------%
% Filename: recover_GDAS_plot_data.m
% Authors: Juan M. Cardenas  
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description:  This script recover the data to plot the figures
%
%-------------------------------------------------------------------------% 

ms1 = 1;    ms2 = 2;    ms3 = 5;    ms4 = 6;    ms = 6;         lw = 1.5;

default_color = [0         0.4470    0.7410 ;
                 0.8500    0.3250    0.0980 ;
                 0.9290    0.6940    0.1250 ;
                 0.4940    0.1840    0.5560 ;
                 0.4660    0.6740    0.1880 ;
                 0.3010    0.7450    0.9330 ;
                 0.6350    0.0780    0.1840 ] ;

markers = {'.','o','x','+','*','s','d','v','^','<','>','p','h',...
    '-*','-o','-s','-^','-v','-p','-h','--*','--o','--s','--^','--v'};
    
%------------------    Save data per dimension    -----------------------%
for p = 1 : length(dim)
    I_Omega = nonzeros(I_Omega_dim(:,1,p));
    for i = 1 : length(N_values)
        if i == 1     
            dV_w(i,p) = length(setdiff(1:K,I_Omega))/length(I_Omega);
            dV_a(i,p) = dV_w(i,p);
            dV_u(i,p) = dV_w(i,p);
        else  
            %--- WLS ---%
            I_lw    = nonzeros(I_ALL_W(:,i-1,Id_Aw(i,p),p));   
            Iw_out1 = setdiff(I_lw,I_Omega);   Iw_out2 = setdiff(I_Omega,I_lw);
            Iw_out  = union(Iw_out1,Iw_out2);       
            dV_w(i,p) = length(Iw_out)/length(I_Omega);
            %--- AWLS ---%
            I_la    = nonzeros(I_ALL_AW(:,i-1,Id_Az(i,p),p));   
            Ia_out1 = setdiff(I_la,I_Omega);   Ia_out2 = setdiff(I_Omega,I_la);
            Ia_out  = union(Ia_out1,Ia_out2);       
            dV_a(i,p) = length(Ia_out)/length(I_Omega);
            %--- LS --- %
            I_lu    = nonzeros(I_ALL_UN(:,i-1,Id_U(i,p),p));   
            Iu_out1 = setdiff(I_lu,I_Omega);   Iu_out2 = setdiff(I_Omega,I_lu);
            Iu_out  = union(Iu_out1,Iu_out2);       
            dV_u(i,p) = length(Iu_out)/length(I_Omega);  
        end        
    end
end

length_dim = length(dim);

Error_matrix = [];
Feval_matrix = [];
Feval_vol_matrix = [];
Rate_matrix = [];
Vol_matrix = [];

for num_dim = 1:length_dim
    
    
    if median_opt == 1
        % Error matrix
        Error_matrix(:,:,num_dim) = [Error_Un_m(:,num_dim), Error_WLS_m(:,num_dim), Error_AWLS_m(:,num_dim), Error_Strat2_m(:,num_dim)];
        
        % C vals and J vals 
        C_vals_matrix(:,:,num_dim) = [C_Un(:,1) C_WLS(:,1) C_AWLS(:,1) C_S2(:,1)];
        J_vals_matrix(:,:,num_dim) = [J_Un(:,1) J_WLS(:,1) J_AWLS(:,1) C_S2(:,1)];
        
        % F eval matrix 
        Feval_matrix(:,:,num_dim) = [Feval_Unif_m(:,num_dim), Feval_WLS_m(:,num_dim), Feval_AWLS_m(:,num_dim), Feval_S2_m(:,num_dim)];
        
        % F eval for Vol matrix 
        Feval_vol_matrix(:,:,num_dim) = [Feval_Unif_m(:,num_dim), Feval_WLS_m(:,num_dim), Feval_AWLS_m(:,num_dim)]; 
        
        % Rejection rate matrix
        Rate_matrix(:,:,num_dim) = [Rate_Unif_m(:,num_dim),  Rate_WLS_m(:,num_dim)];
        
        
        
    else
        % Error matrix
        Error_matrix(:,:,num_dim) = [Error_Un_mn(:,num_dim), Error_WLS_mn(:,num_dim), Error_AWLS_mn(:,num_dim), Error_Strat2_mn(:,num_dim)];
        
        % C vals 
        C_vals_matrix(:,:,num_dim) = [C_Un_mn(:,1) C_WLS_mn(:,1) C_AWLS_mn(:,1) C_S2_mn(:,1)];        
        J_vals_matrix(:,:,num_dim) = [J_Un_mn(:,1) J_WLS_mn(:,1) J_AWLS_mn(:,1) C_S2_mn(:,1)];         
        
        % F eval matrix
        Feval_matrix(:,:,num_dim) = [Feval_Unif_mn(:,num_dim), Feval_WLS_mn(:,num_dim), Feval_AWLS_mn(:,num_dim), Feval_S2_mn(:,num_dim)]; 
        
        % F eval for Vol matrix
        Feval_vol_matrix(:,:,num_dim) = [Feval_Unif_mn(:,num_dim), Feval_WLS_mn(:,num_dim), Feval_AWLS_mn(:,num_dim)];
        
        % Rejectio rate matrix
        Rate_matrix(:,:,num_dim) = [Rate_Unif_mn(:,num_dim),  Rate_WLS_mn(:,num_dim)];
        
    end
    
    % Volume matrix 
    Vol_matrix(:,:,num_dim)  = [dV_u(:,num_dim), dV_w(:,num_dim), dV_a(:,num_dim)]*100; 
    
    % N values matrix 
    N_values_matrix(:,1:4,num_dim) = N_values(:,num_dim).*ones(length(N_values(:,1)),4);
    
    % M values matrix 
    M_values_matrix(:,1:4,num_dim) = M_values(:,num_dim).*ones(length(M_values(:,1)),4);
    
    % numer of iterations
    l_steps = [1:length(N_values(:,1))]'.*ones(length(M_values(:,1)),4);
  
end
  