%-------------------------------------------------------------------------%
% Filename: Domain_examples.m
% Authors: Juan M. Cardenas  
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description:  This script computes the volume for a domain of interest 
% for a function in several dimensions. Also, 
%
% Inputs:
% f - function
% K - size grid
% dim - vector of dimensions
% domain - type of domain 
% tol - parameter for domain 1: f>= tol
% a,b - parameter for domain 2: a<= f <= b
% 
% Output: 
% K_Omega - size of grid each learned domain for each dimension
% Vol_Omega - Volume for the domain of interest 
% plot - plot of volume through dimensions
%-------------------------------------------------------------------------% 
close all; clear all ; clc;

f      = 13;
K      = 30000;
dim    = [2 3 4 5 6 7 8 10 11 12 13 14 15]; 
domain = 2;
tol    = 0;
a      = 1.8/10;                   % Domain 2 parameter
b      = 0.72;                     % Domain 2 parameter 
space  = '';
indices_D = [1:K]';
%-------------------------------------------------------------------------%
% Plot parameters 
ms1 = 1;    ms2 = 2;    lw = 1;
default_color = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.9290    0.6940    0.1250 ;
    0.4940    0.1840    0.5560 ;
    0.4660    0.6740    0.1880 ;
    0.3010    0.7450    0.9330 ;
    0.6350    0.0780    0.1840 ] ;

markers = {'.','o','x','+','*','s','d','v','^','<','>','p','h',...
    '-*','-o','-s','-^','-v','-+','-x','--*','--o','--s','--^','--v'};
%-------------------------------------------------------------------------%
Vol_Omega = [];
K_Omega   = [];

for p = 1:length(dim)
    
    d = dim(p);
    
    % Grid 
    Z = 2*rand(K,d)-1;

    % Plot Grid Z
    if d == 2
        fig1 = figure();
        plot(Z(:,1),Z(:,2),markers{1},'color',default_color(7,:),...
            'MarkerSize',ms1);
        T = title('Domain','interpreter','latex');
    elseif d == 3 
        fig1 = figure();
        plot3(Z(:,1),Z(:,2),Z(:,3),'.')
        T = title('Domain','interpreter','latex');
    else 
        disp('Dimension is bigger than 3')
    end
    
    % Function F   
    switch f
        case 1
            % high-dim: Annulus
            Y      = sum(Z.^2,2);
            f_grid =  log(8*Y)-(2)*Y;
        case 2          %  rare, old function 8
            y1 = Z(:,1);    y2 = Z(:,2);
            f_grid = log(8*(y1.^2 + y2.^2))-3*y1.^2 +y2.^2 +cos(2*y1).*sin(y2)-exp(-20*((y2-0.5).^2));    
        case 3          % d-dimensions Prism
            y1 = Z(:,1);    y2 = Z(:,2);
            f_grid =  (4 - 1./(y1.^2 + y2.^2)).*exp(-(sum(Z,2)/d));         
        case 4
            Y = Z;
            for z = 1 : K
                for k = (d/2+1) : d
                    f_1(k-d/2) = cos(16*Y(z,k)/2^k);
                end
                f_1 = prod(f_1);
                for j = 1 : (d/2)
                    f_2(j) = 1 - (Y(z,j)/2^j); 
                end
                f_2 = prod(f_2);
                f_grid(z,1) = f_1*f_2*(log(8*sum(Y(z,:).^2,2))-2*(sum(Y(z,:).^2,2)));
            end
        case 5
            Y = Z;
            for z = 1 : K
                for k = (d/2+1) : d
                    f_1(k-d/2) = cos(16*Y(z,k)/(2^k));
                end
                f_1 = prod(f_1);
                for j = 1 : (d/2)
                    f_2(j) = 1 - (Y(z,j)/4^j); 
                end
                f_2 = prod(f_2);
                f_grid(z) = f_1/f_2;
            end  
        case 6
            Y = sum(Z,2);
            f_grid = exp(-Y/(2*d)); 
        case 7 
            Y_1 = sum(Z,2);         Y_2 = prod(Z,2);
            f_grid = exp(-Y_1/(2*d)).*sin(-pi*Y_1); 
        case 8
            Y_1 = sum(Z,2);         Y_2 = prod(Z,2);
            f_grid = exp(-Y_2/(d)).*sin(2*pi*Y_1).*log(Y_1/4);
        case 9
            Y_1 = sum(Z,2);         Y_2 = prod(Z,2);
            f_grid = exp(-Y_2/d)./(Y_1+1/d) - Y_2;     
        case 10 
            f_grid = fun_one(Z);
        case 11 
            Y = Z;
            for z = 1:K
                for k = 1:ceil(d/2)
                    f_1(k) = 1 + (4^k)*(Y(z,k)^2);
                end
                f_1 = prod(f_1);
                for j = ceil(d/2)+1 : d
                    f_2(j-ceil(d/2)) = 100+5*Y(z,j);
                end
                f_2 = prod(f_2);
                f_grid(z,1) = (f_1/f_2)^(1/d);
            end
        case 12
            % High-dim: New Annulus
            Y      = sum(Z.^2,2);
            f_grid = log((16/d)*Y)-(4/d)*Y;
        case 13   
            ind    = [1:d] + 1;
            v_1    = (Z + (((-1).^ind)./ind)).^2;
            f_grid = prod( (d/4)./((d/4) + v_1),2);            
    end

    % Compute Omega
    if domain == 1
        I_Omega = indices_D(f_grid >= tol);
    else
        % domain 2 
        %I_Omega = indices_D(a <= f_grid <= b);
        I_om1 = indices_D(a<= f_grid);
        I_om2 = indices_D(f_grid <= b);
        I_Omega = intersect(I_om1,I_om2);
    end
    
    % Plot Omega 
    if d == 2
        fig2 = figure();
        plot(Z(I_Omega,1),Z(I_Omega,2),markers{1},'color',...
            default_color(4,:),'MarkerSize',ms1);
        T = title('$\Omega$','interpreter','latex');
    elseif d == 3
        fig2 = figure();
        plot3(Z(I_Omega,1),Z(I_Omega,2),Z(I_Omega,3),markers{1},...
            'color',default_color(4,:),'MarkerSize',ms1);
        T = title('$\Omega$','interpreter','latex');
    else 
        disp('The dimension is bigger than 3')
    end

    ax          = gca;
    ax.FontSize = 12;                       ax.LineWidth = 1;
    ax.XTick    = [-1 -0.5 0 0.5 1];        xlim([-1 1]);
    ax.YTick    = [-1 -0.5 0 0.5 1];        ylim([-1 1]);

    K_Omega(p)   = length(I_Omega);
    Vol_Omega(p) = (K_Omega(p)/K)*100;
end

%-------------------------------------------------------------------------%
Show_results = [K_Omega', Vol_Omega']
figure()
plot(dim,Vol_Omega,'r-o',dim,70*ones(length(Vol_Omega)),'b--*')
xlabel('dimension','interpreter','latex')
ylabel('Volume of $\Omega$','interpreter','latex')
grid on 