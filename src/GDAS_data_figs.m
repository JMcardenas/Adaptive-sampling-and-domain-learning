%-------------------------------------------------------------------------%
% Filename: GDAS_data_figs.m
% Authors: Juan M. Cardenas and Ben Adcock
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description: This script computes ASUD for weighted-LS, augmented 
% weighted-LS, Monte Carlo uniform sampling LS, and ASGD weighted-LS
%
% Inputs: 
% index_type - type of index set (1-hyperbolic,0-total degree)
% dim - vector of dimensions
% r_f - ratio for functions, in particular f=f_1
% Trials - number of trials
% save_file - parameter to save files (1-true,0-false)
% N_max - maximum size of index set
% f - type of function (from 1 to 4)
% domain - type of domain: 1- f>= tol, 2- a<= f <= b
% tol - parameter for domain of interest: f >= tol
% a,b - parameters for domain of interest: a<= f <= b
%
% Outputs: 
% Error - Relative error 
% s_min - minimum singular value 
% C vals - Constant C values 
% I_all - Indices of Grid Z_i: I_all    
% Rate_f - rate of functions evaluations
% Feval - functions evaluations
% dV - Volume of learned domain
% J vals - Constant J values 
%-------------------------------------------------------------------------%
%% Set up
tic
space=' ';                                             

% Radio for f=3: Normal setup r_1=1/2(80% Uniform)| r_2=7/10(60%)|r_3=7/8(40%)  
if r_f  == 1		  
	r_d = 1/2;
elseif r_f == 2
	r_d = 7/10;
else  	
	r_d = 7/8;
end
 
n_values  = [];
K         = 30*N_max;        % size of the grid
indices_D = [1:K]';          % Indices of D domain

for p = 1 : length(dim)
    
    d = dim(p);
    %% Pre-processing
    %---------------------------------------------------------------------%
    %---                 Generate Grid Z over Omega                    ---%
    %---------------------------------------------------------------------%
    % Unitary Hypercube 

    Z = 2*rand(K,d)-1;
    
    if plt_domain == 1
        if d == 2
            plot(Z(:,1),Z(:,2),'r.');
        elseif d == 3
            plot3(Z(:,1),Z(:,2),Z(:,3),'r.');
        else
            disp('dimension bigger than 3');
        end
    end    
    
    %---------------------------------------------------------------------%
    %---                        Index set                              ---%
    %---------------------------------------------------------------------%  
    HCsize = 0;         n_0 = 0; 
    
    while HCsize < N_max
        n_0 = n_0 + 1;
        HCsize = length(HC_index(n_0,d));
    end
    nmax = n_0 -1;      Nmax = length(HC_index(nmax,d));
    
    if d < 15
        n_values_p(:,p) = round(linspace(1,nmax,max_div));
        n_values        = n_values_p(2:end,p);                  % Avoid n=1
    else  
        n_values = round(linspace(1,nmax,max_div-1));
    end
    
    % Re-order index set      
    for i = 1 : length(n_values)    
        I = HC_index(n_values(i),d);
        N_values(i,p) = length(I);
        if i == 1
            J = I;
        else
            C = setdiff(I',J','rows');              % Index in I\J
            C = C';                                   
            J = [J  C];                             % Adding C to J indices
        end
    end
    
    %---------------------------------------------------------------------%
    %---                Generate Measurement matrix                    ---% 
    %---------------------------------------------------------------------%   
    B = zeros(K,Nmax);              
    
    for i = 1:K                     
        z = Z(i,:);                     
        L = LegMat(z',nmax+1);      
        for j = 1:Nmax             
            % Loop over the index set
            Lij = zeros(d,1); 
            for k = 1:d                   
                % Loop over the components 
                Lij(k,1) = L(k,J(k,j)+1); 
            end
            % Tensor product-type
            B(i,j) = prod(Lij);           
        end
    end          
    
    %---------------------------------------------------------------------%    
    %---                Function right-hand side                       ---%
    %---------------------------------------------------------------------%    

    switch f
        case 1          
            % high-dim: Prism
            y1     = Z(:,1);    
            y2     = Z(:,2);
            f_grid = ((1/r_d)^2 - 1./(y1.^2 + y2.^2)).*exp(-(sum(Z,2)/(2*d)));        
        case 2  
            % High-dim: Annulus
            Y      = sum(Z.^2,2);
            f_grid = log(8*Y)-(2)*Y;             
        case 3          
            % High-dim: New f_1
            Y      = sum(Z.^2,2);
            f_grid = g_dim(d)*log((16/d)*Y)-(4/d)*Y;              
        case 4   
            ind    = [1:d] + 1;
            v_1    = (Z + (((-1).^ind)./ind)).^2;
            f_grid = prod( (d/4)./((d/4) + v_1),2);            
    end
    
    %---------------------------------------------------------------------%
    %---                        Find Omega                             ---%
    %---------------------------------------------------------------------%     

    if domain == 1
        I_Omega = indices_D(f_grid >= tol);
    else
        % domain 2 
        %I_Omega = indices_D(a <= f_grid <= b);
        I_om1 = indices_D(a<= f_grid);
        I_om2 = indices_D(f_grid <= b);
        I_Omega = intersect(I_om1,I_om2);
    end
    K_Omega = length(I_Omega);

    % Plot Omega
    if plt_domain == 1
        if d == 2
            plot(Z(I_Omega,1),Z(I_Omega,2),'b.');
        elseif d == 3
            plot3(Z(I_Omega,1),Z(I_Omega,2),Z(I_Omega,3),'b.');
        else
            disp('dimension bigger than 3');
        end
    end 

    disp(['------------------------- START --------------------------------']);
    disp([ 'Size Grid: ',num2str(K),space,'| Pts Omega: ',num2str(K_Omega),...
            '| Domain: ',num2str(domain),space,'| function: ',num2str(f)])
    disp(['----------------------------------------------------------------']);
           
    %---------------------------------------------------------------------%    
    %---                          M points                             ---%
    %---------------------------------------------------------------------%
    k_pts         = round(log(N_values(:,p)));      % M = kN 
    M_values(:,p) = k_pts.*N_values(:,p);
    M_val         = [0; M_values(:,p)];              
    N_val         = [0; N_values(:,p)];
    
    %% Loop over Trials
    I_ALL = zeros(K,max_div -1,Trials);

    % ASGD-WLS: Orthogonalize over P_l
    B_s2        = (1/sqrt(K_Omega))*B(I_Omega,:);
    [Q_s2,R_s2] = qr(B_s2,0);                      % QR for ASGD-WLS
    mu_s2       = abs(Q_s2.^2);                    % sampling measure

     for t = 1 : Trials
        % Auxiliary variables
        s_min_Aw = [];  Error_w  = [];  Re_w = [];  I_w = [];  
        s_min_Aa = [];  Error_a  = [];  R_a  = [];  I_a = [];  I_in_Om = [];
        s_min_s2 = [];  Error_s2 = [];              I_s = [];     
        s_min_un = [];  Error_un = [];  R_u = [];   I_u = [];        
   
        F_ev_w    = zeros(length(N_values(:,p)),1);           
        F_ev_a    = F_ev_w;              
        F_ev_un   = F_ev_w;           

        for l = 1 : length(N_values(:,p))
            %%  ASUD - WLS

            N = N_values(l,p);
            M = M_values(l,p);

            % Orthogonalize P_l over Z_l
            if l == 1
                B_w  = B(:,1:N);
                I_iw = 1:K ;  
            else 
                B_w  = B(I_jw,1:N);
                I_iw = I_jw;
            end

            K_iw      = length(I_iw);    
            B_w       = sqrt(1/K_iw)*B_w;        % Montecarlo factor
            [Q_w,R_w] = qr(B_w,0);               % Reduce QR-decom  

            % Sampling measure mu 
            mu_w = abs(Q_w.^2);  
            
            % Draw pts from psi_l
            for j = 1:N
                if j <= N_val(l)
                    k = k_pts(l)-k_pts(l-1);
                else
                    k = k_pts(l);
                end
                i = 0;
                while i < k                
                    z          = datasample(I_iw,1,'Replace',true,'Weights',mu_w(:,j));   
                    F_ev_w(l) = F_ev_w(l) + 1;  
                    if const_dom(f_grid(z),domain,tol,a,b) == 1
                        I_w = [I_w; z];     
                        i    = i + 1;
                    else
                        Re_w = [Re_w; z]; 
                    end
                end
            end 
         
            % Find M_index in I set
            I_M = [];
            for i_m = 1 : M                             
                I_M = [I_M ; find(I_iw == I_w(i_m))];   
            end
            
            % Compute approximation
            W   = ((1/N)*(sum(mu_w(I_M,1:N)'))).^(-1);  % Compute the weights 
            W_w = diag(sqrt(W));                        
            A_w = (1/sqrt(M))*W_w*Q_w(I_M,1:N); 
            b_w = (1/sqrt(K_iw*M))*W_w*f_grid(I_w); 
               
            % Compute LS-approx
            c_w         = A_w\b_w; 
            s_min_Aw(l) = min(svd(A_w));     
            
            % J constant
            J_const_indices = [];
            for i = 1 : M
                J_const_indices = [J_const_indices, find(I_Omega == I_w(i))]; 
            end
            A_tilde = (1/sqrt(M))*W_w*Q_s2(J_const_indices,1:N);
            J_Aw_temp(l) = 1/min(svd(A_tilde));
            
            C(1:length(c_w),l) = c_w; 

            % Relative error    
            f_tilde_w  = B(:,1:N)*(R_w\c_w);    
            Error_w(l) = norm(f_grid(I_Omega)-f_tilde_w(I_Omega))/norm(f_grid(I_Omega));
           
            %-------------------------------------------------------------%
            %---                 Compute Next Grid                     ---%
            %-------------------------------------------------------------%
            if domain == 1
               I_z = indices_D(f_tilde_w >= tol);
            else
               % domain 2
               %I_z = indices_D(a <= f_tilde_w <= b);
                I_z1 = indices_D(a <= f_tilde_w);
                I_z2 = indices_D(f_tilde_w <= b);
                I_z = intersect(I_z1,I_z2);
            end            
            I_jw = setdiff(I_z,Re_w);                           

            % Adding M pts to Z_l
            I_jw = union(I_w,I_jw);       

            % Save Grid 
            I_all_w(1:length(I_jw),l) = I_jw;     
            
            
            %-------------------------------------------------------------%
            % print info 
            disp(['d = ',num2str(d),space,'trial = ',num2str(t),space,'iter = ',num2str(l),space,...
                  'num pts = ',num2str(M),space,'num basis = ',num2str(N),space,'F eval = ',...
                  num2str(F_ev_w(l)),space,'Error = ',num2str(Error_w(l),'%10.2e'),space,...
                  'C = ',num2str(1/s_min_Aw(l),'%10.2e'),space,'J = ',num2str(J_Aw_temp(l),'%10.2e'),...
                  space,'Method = ASUD-LS']);
            
            %%  ASUD - Augmented WLS

            % Orthogonalize P_i over Z_i 
            if l == 1
                B_a  = B(:,1:N);
                I_ia = [1:K]';  
            else 
                B_a  = B(I_ja,1:N);
                I_ia = I_ja;
            end
            K_ia      = length(I_ia);          
            B_a       = sqrt(1/K_ia)*B_a;        % Montecarlo factor
            [Q_a,R_a] = qr(B_a,0);               % Reduce QR-decom  
            
            % Sampling measure
            mu_a = abs(Q_a.^2);    
    
            % Draw pts from psi_l
            for j = 1:N
                if j <= N_val(l)
                    k = k_pts(l)-k_pts(l-1);
                else
                    k = k_pts(l);
                end  
                I_ad = datasample(I_ia,k,'Replace',true,'Weights',mu_a(:,j));   
                I_a  = [I_a; I_ad];     
            end
 
            % Find M_index in I set
            I_M = [];
            for i_m = 1 : M                           
                I_M = [I_M ; find(I_ia == I_a(i_m))];   
            end
            
            % Compute approximation
            W   = ((1/N)*(sum(mu_a(I_M,1:N)'))).^(-1);  % Compute the weights 
            W_a = diag(sqrt(W));                        
            A_a = (1/sqrt(M))*W_a*Q_a(I_M,1:N); 
            b_a = (1/sqrt(K_ia*M))*W_a*f_grid(I_a);

            % Compute LS-approx
            c_a                = A_a\b_a; 
            s_min_Aa(l)        = min(svd(A_a));           
            C(1:length(c_a),l) = c_a;
            
            % J constant
            J_const_indices = [];
            i_indices = [];
            I_inside_Omega = [];
            for i = 1 : M
                J_const_indices = [J_const_indices; find(I_Omega == I_a(i))]; 
                if length(find(I_Omega == I_a(i))) == 1
                    i_indices = [i_indices; i];
                end
            end
            
            for i = 1 : length(i_indices)
                I_inside_Omega = [I_inside_Omega; find(I_ia == I_a(i_indices(i)))];
            end
            W_inside_Omega = ((1/N)*(sum(mu_a(I_inside_Omega,1:N)'))).^(-1);  % Compute the weights 
            W_inside_Omega = diag(sqrt(W_inside_Omega)); 
            A_tilde = (1/sqrt(length(i_indices)))*W_inside_Omega*Q_s2(J_const_indices,1:N);
            J_Aa_temp(l) = 1/min(svd(A_tilde));
            
            % Relative error    
            f_tilde_a  = B(:,1:N)*(R_a\c_a);    
            Error_a(l) = norm(f_grid(I_Omega)-f_tilde_a(I_Omega))/norm(f_grid(I_Omega));

            %-------------------------------------------------------------%
            %---                 Compute Next Grid                     ---%
            %-------------------------------------------------------------%
            if domain == 1
                I_z = indices_D(f_tilde_a >= tol);
            else
                % domain 2
                %I_z = indices_D(a <= f_tilde_a <= b);
                I_z1 = indices_D(a <= f_tilde_a);
                I_z2 = indices_D(f_tilde_a <= b);
                I_z = intersect(I_z1,I_z2);
            end
            
            % Add M pts in Omega 
            if domain == 1
                I_in_Om = I_a( (f_grid(I_a)>= 0) );
                Re_a    = I_a( (f_grid(I_a)< 0) );
            else 
                % domain 2 
                I_in_Om1 = I_a(a<= f_grid(I_a));
                I_in_Om2 = I_a(f_grid(I_a) <= b);
                I_in_Om = intersect(I_in_Om1,I_in_Om2);
                Re_a    = setdiff(I_a,I_in_Om);                
                %Re_a    = I_a( (f_grid<a && b<f_grid) );
                %I_in_Om = I_a(  (a <= f_grid(I_a) <= b) );
            end
            count_I_in_omega(l) = length(I_in_Om);
    
            % Adding M pts to Z_l
            I_ja = union(I_a,I_z);     

            % Taking off rejected points 
            I_ja_final = setdiff(I_ja,Re_a);

            % Save Grid_I
            I_all_a(1:length(I_ja_final),l) = I_ja_final;  
            J_exc_all(1:length(Re_a),l)     = Re_a;
            
            %-------------------------------------------------------------%
            % print info 
            disp(['d = ',num2str(d),space,'trial = ',num2str(t),space,'iter = ',num2str(l),space,...
                  'num pts = ',num2str(M),space,'num basis = ',num2str(N),space,'F eval = ',...
                  num2str(M),space,'Error = ',num2str(Error_a(l),'%10.2e'),space,'C = ',...
                  num2str(1/s_min_Aa(l),'%10.2e'),space,'J = ',num2str(J_Aa_temp(l),'%10.2e'),...
                  space,'Method = ASUD-ALS']);         
  
            %%  ASGD - WLS  
    
            % Draw pts from psi_l
            for j = 1:N
                if j <= N_val(l)
                    k = k_pts(l)-k_pts(l-1);
                else
                    k = k_pts(l);
                end  
                I_ad = datasample([1:K_Omega]',k,'Replace',true,'Weights',mu_s2(:,j));   
                I_s  = [I_s; I_ad];     
            end
    
            % Find Is points in I_Omega 
            I_s2 = I_Omega(I_s);

            % Compute LS-approx
            Ws   = ((1/N)*(sum(mu_s2(I_s,1:N)'))).^(-1);   % Compute the weights  
            W_Ds = diag(sqrt(Ws));                         
            As   = (1/sqrt(M))*W_Ds*Q_s2(I_s,1:N);            
            bs   = (1/sqrt(M*K_Omega))*W_Ds*f_grid(I_s2);

            c_s2        = As\bs;
            s_min_s2(l) = min(svd(As));

            % Relative error
            f_tilde_s2  = sqrt(K_Omega)*Q_s2(:,1:N)*c_s2;
            Error_s2(l) = norm(f_grid(I_Omega) - f_tilde_s2)/norm(f_grid(I_Omega));

            %-------------------------------------------------------------%
            % print info 
            disp(['d = ',num2str(d),space,'trial = ',num2str(t),space,'iter = ',num2str(l),space,...
                  'num pts = ',num2str(M),space,'num basis = ',num2str(N),space,'F eval = ',...
                  num2str(M),space,'Error = ',num2str(Error_s2(l),'%10.2e'),space,'C = ',...
                  num2str(1/s_min_s2(l),'%10.2e'),space,'J = ',num2str(1/s_min_s2(l),'%10.2e'),...
                  space,'Method = ASGD-LS']);           

            %% Uniform - LS
    
            % Draw M_l pts uniform 
            while length(I_u) < M 
                z          = datasample(1:K,1,'Replace',true);  
                F_ev_un(l) = F_ev_un(l) + 1;
                if const_dom(f_grid(z),domain,tol,a,b) == 1
                    I_u = [I_u; z];     
                else
                    R_u = [R_u; z];
                end 
            end  

            % Compute system 
            A_u = (1/sqrt(M))*B(I_u,1:N);
            b_u = (1/sqrt(M))*f_grid(I_u);

            % Compute LS-Approx
            c_u         = A_u\b_u; 
            s_min_un(l) = min(svd(A_u));                  
    
            C_u(1:length(c_u),l) = c_u; 
        
            % J constant
            J_const_indices = [];
            for i = 1 : M
                J_const_indices = [J_const_indices, find(I_Omega == I_u(i))]; 
            end
            A_tilde = (1/sqrt(M))*Q_s2(J_const_indices,1:N);
            J_Au_temp(l) = 1/min(svd(A_tilde));

            % Compute error
            f_tilde_un  = B(:,1:N)*c_u;                    
            Error_un(l) = norm(f_grid(I_Omega) - f_tilde_un(I_Omega))/norm(f_grid(I_Omega));

            %-------------------------------------------------------------%
            %---                 Compute Next Grid                     ---%
            %-------------------------------------------------------------%
            if domain == 1
                I_unz = indices_D(f_tilde_un >= tol);
            else
                % domain 2
                %I_unz = indices_D(a <= f_tilde_un <= b);
                I_un1 = indices_D(a <= f_tilde_un);
                I_un2 = indices_D(f_tilde_un <= b);
                I_unz = intersect(I_un1,I_un2);
            end
            % Extract rejected pts
            I_unj = setdiff(I_unz,R_u);
            
            % Adding M pts to Z_l
            I_unj = union(I_u,I_unj);  
            
            % Save Grid 
            I_all_un(1:length(I_unj),l) = I_unj;    
            
           %-------------------------------------------------------------%
            % print info 
            disp(['d = ',num2str(d),space,'trial = ',num2str(t),space,'iter = ',num2str(l),space,...
                  'num pts = ',num2str(M),space,'num basis = ',num2str(N),space,'F eval = ',...
                  num2str(F_ev_un(l)),space,'Error = ',num2str(Error_un(l),'%10.2e'),space,'C = ',...
                  num2str(1/s_min_un(l),'%10.2e'),space,'J = ',num2str(J_Au_temp(l),'%10.2e'),...
                  space,'Method = MC-LS']);
              
        end            

        %-----------------------------------------------------------------%
        %--------                 Save Info                     ----------%
        %-----------------------------------------------------------------%
        % ASUD-WLS
        Error_WLS(:,t,p) = Error_w;            s_min_WLS(:,t,p) = s_min_Aw;    
        Feval_WLS(:,t,p) = cumsum(F_ev_w);     J_const_Aw(:,t,p) = J_Aw_temp;

        I_ALL_W(1:length(I_all_w),:,t,p)  = I_all_w;
        I_WLS_M(1:max(M_values(:,p)),t,p) = I_w; 

        % ASUD-AWLS
        Error_AWLS(:,t,p) = Error_a;           s_min_AWLS(:,t,p) = s_min_Aa;     
                                               J_const_Aa(:,t,p) = J_Aa_temp;
                                               
        I_ALL_AW(1:length(I_all_a),:,t,p)   = I_all_a;
        I_ALL_AR(1:length(J_exc_all),:,t,p) = J_exc_all;
        I_WLS_AM(1:length(I_in_Om),t,p)     = I_in_Om; 
        Count_I_m(1:l,t,p)                  = count_I_in_omega; 

        % ASGD-WLS
        Error_Strat2(:,t,p) = Error_s2;       s_min_strat2(:,t,p) = s_min_s2;    

        % MC-LS 
        Error_Unif(:,t,p) = Error_un;         s_min_Unif(:,t,p) = s_min_un;    
        Feval_Unif(:,t,p) = cumsum(F_ev_un);  J_const_Au(:,t,p) = J_Au_temp;

        I_ALL_UN(1:length(I_all_un),:,t,p) = I_all_un;
        I_Unif_M(1:max(M_values(:,p)),t,p) = I_u;

        % Save grid
        I_Omega_dim(1:length(I_Omega),t,p)  = I_Omega;  
        Z_l_dim(:,1:d,p)                    = Z;    
    end
end

t1 = toc;

for p = 1:length(dim)
    % Median of data 
    Error_WLS_m(:,p)     = median(Error_WLS(:,:,p),2);
    Error_AWLS_m(:,p)    = median(Error_AWLS(:,:,p),2);
    Error_Strat2_m(:,p)  = median(Error_Strat2(:,:,p),2);
    Error_Un_m(:,p)      = median(Error_Unif(:,:,p),2);

    Feval_WLS_m(:,p)     = median(Feval_WLS(:,:,p),2);
    Feval_AWLS_m(:,p)    = M_values(:,p);
    Feval_S2_m(:,p)      = M_values(:,p);
    Feval_Unif_m(:,p)    = median(Feval_Unif(:,:,p),2);


    sv_min_WLS(:,p)      = median(s_min_WLS(:,:,p),2);
    sv_min_AWLS(:,p)     = median(s_min_AWLS(:,:,p),2); 
    sv_min_S2(:,p)       = median(s_min_strat2(:,:,p),2);     
    sv_min_Unif(:,p)     = median(s_min_Unif(:,:,p),2);     

    C_WLS(:,p)           = 1./sv_min_WLS(:,p);
    C_AWLS(:,p)          = 1./sv_min_AWLS(:,p); 
    C_S2(:,p)            = 1./sv_min_S2(:,p);
    C_Un(:,p)            = 1./sv_min_Unif(:,p);
    
    J_WLS(:,p)           = median(J_const_Aw(:,:,p),2);
    J_AWLS(:,p)          = median(J_const_Aa(:,:,p),2);
    J_Un(:,p)            = median(J_const_Au(:,:,p),2);
    
    % Mean of data 
    Error_WLS_mn(:,p)    = mean(Error_WLS(:,:,p),2);
    Error_AWLS_mn(:,p)   = mean(Error_AWLS(:,:,p),2);
    Error_Strat2_mn(:,p) = mean(Error_Strat2(:,:,p),2);
    Error_Un_mn(:,p)     = mean(Error_Unif(:,:,p),2);


    Feval_WLS_mn(:,p)    = mean(Feval_WLS(:,:,p),2);
    Feval_AWLS_mn(:,p)   = M_values(:,p);
    Feval_S2_mn(:,p)     = M_values(:,p);
    Feval_Unif_mn(:,p)   = mean(Feval_Unif(:,:,p),2);


    sv_min_WLS_mn(:,p)   = mean(s_min_WLS(:,:,p),2);
    sv_min_AWLS_mn(:,p)  = mean(s_min_AWLS(:,:,p),2); 
    sv_min_S2_mn(:,p)    = mean(s_min_strat2(:,:,p),2);     
    sv_min_Unif_mn(:,p)  = mean(s_min_Unif(:,:,p),2);     

    C_WLS_mn(:,p)        = 1./sv_min_WLS_mn(:,p);
    C_AWLS_mn(:,p)       = 1./sv_min_AWLS_mn(:,p); 
    C_S2_mn(:,p)         = 1./sv_min_S2_mn(:,p);
    C_Un_mn(:,p)         = 1./sv_min_Unif_mn(:,p);    

    J_WLS_mn(:,p)        = mean(J_const_Aw(:,:,p),2);
    J_AWLS_mn(:,p)       = mean(J_const_Aa(:,:,p),2);
    J_Un_mn(:,p)         = mean(J_const_Au(:,:,p),2);    
end

%--- Find the median's component ---% 
Vec = ones(Trials ,1);
for p = 1 : length(dim)
    for j = 1: max_div - 1
        Idx        = knnsearch(s_min_WLS(j,:,p)',sv_min_WLS(j,p)*Vec);
        Id_Aw(j,p) = Idx(1);
        Idy        = knnsearch(s_min_Unif(j,:,p)',sv_min_Unif(j,p)*Vec);
        Id_U(j,p)  = Idy(1);
        Idz        = knnsearch(s_min_AWLS(j,:,p)',sv_min_AWLS(j,p)*Vec);
        Id_Az(j,p) = Idz(1);      
    end
end

for p = 1 : length(dim)
    I_Omega = nonzeros(I_Omega_dim(:,1,p));
    for i = 1 : length(N_values(:,p))
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

%---  Rate F-evaluations  ---%  
 
% Median
Rate_WLS_m   = (Feval_WLS_m  - M_values)./Feval_WLS_m; 
Rate_Unif_m  = (Feval_Unif_m - M_values)./Feval_Unif_m;

% Mean
Rate_WLS_mn  = (Feval_WLS_mn  - M_values)./Feval_WLS_mn; 
Rate_Unif_mn = (Feval_Unif_mn - M_values)./Feval_Unif_mn; 

%% Save info 

if save_file == 1 
    if dim == 2     
        if fig_num == 7 
            namefile = ['../data/fig' ' ' sprintf('%d',fig_num) '/' sprintf('Data_GDAS_fig_%d_f_%d_dim2.mat',fig_num,f)];
        else
            namefile = ['../data/fig' ' ' sprintf('%d',fig_num) '/' sprintf('Data_GDAS_fig_%d_dim2.mat',fig_num)];
        end
    else 
       namefile = ['../data/fig' ' ' sprintf('%d',fig_num) '/'  sprintf('Data_GDAS_fig_%d.mat',fig_num)];
    end
    save(namefile,'N_values','M_values',...
    'Error_WLS','s_min_WLS','Feval_WLS','Error_WLS_m','Feval_WLS_m',...
    'sv_min_WLS','C_WLS','Error_WLS_mn','Feval_WLS_mn','sv_min_WLS_mn',...
    'C_WLS_mn',...
    'Error_AWLS','s_min_AWLS','Error_AWLS_m','Feval_AWLS_m',...
    'sv_min_AWLS','C_AWLS','Error_AWLS_mn','Feval_AWLS_mn',...
    'sv_min_AWLS_mn','C_AWLS_mn',...
    'Error_Unif','s_min_Unif','Feval_Unif','Error_Un_m','Feval_Unif_m',...
    'sv_min_Unif','C_Un','Error_Un_mn','Feval_Unif_mn','sv_min_Unif_mn',...
    'C_Un_mn',...
    'Error_s2','s_min_strat2','Error_Strat2_m','Feval_S2_m','sv_min_S2',...
    'C_S2','Error_Strat2_mn','Feval_S2_mn','sv_min_S2_mn','C_S2_mn',...
    'Rate_WLS_m','Rate_Unif_m','Rate_WLS_mn','Rate_Unif_mn',...
    'I_ALL_W','I_ALL_AW','I_WLS_M','I_WLS_AM','I_ALL_UN','I_Unif_M',...
    'I_Omega','Z','Count_I_m',...
    'I_Omega_dim','dV_w','dV_a','dV_u','Id_Aw','Id_Az','Id_U','Z_l_dim',...
    'K','Trials','max_div','N_max','tol','f','r_f','dim','t1',...
    'J_const_Aw','J_const_Aa','J_const_Au','J_WLS','J_AWLS','J_Un',...
    'J_WLS_mn','J_AWLS_mn','J_Un_mn');
end
