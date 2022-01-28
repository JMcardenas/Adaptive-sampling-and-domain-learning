% set fonts run_GDAS_plot_figs_3_to_6

            if col_num == 1
                set(gca,'xscale','log');
                set(gca,'yscale','log');
                h = legend('MC-LS','ASUD-LS','ASUD-ALS','ASGD-LS');
                
            elseif col_num == 2
                set(gca,'xscale','log'); 
                h = legend('MC-LS','ASUD-LS','ASUD-ALS');
                
            elseif col_num == 3
                h = legend('MC-LS','ASUD-LS');              
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