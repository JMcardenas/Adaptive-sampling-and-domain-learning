%-------------------------------------------------------------------------%
% Filename: beautify_plot.m
% Authors: Juan M. Cardenas  
% Part of the paper "An adaptive sampling and domain learning strategy for  
% multivariate function approximation on unknown domains"
%
% Description:  This script modify the plot for position purpose 
%
% Inputs:
% gcf - figure 
% 
% Output: 
% figure 
%-------------------------------------------------------------------------% 
gcf;

set(gcf,'Color',[1 1 1]);

font_size = 20;

set(findall(gcf,'type','text'),'fontSize',font_size);
set(findall(gcf,'type','axes'),'fontsize',font_size);
set(gcf,'DefaultTextFontSize',font_size);

%Position plot at left hand corner with width 10 and height 10.
set(gcf, 'PaperPosition', [0 0 10 10]); 

%Set the paper to have width 10 and height 10.
set(gcf, 'PaperSize', [10 10]);
set(gcf,'Position',[1 1 800 800]);

ax = gca;
ax.Position = [0.1 0.1 0.83 0.83];

nplots = size(ax.Children,1);
 
xmin = 1e16;
xmax = 1e-16;
ymin = 1e16;
ymax = 1e-16;

for j = 1:nplots
%     ax.Children(j).DisplayName = ['data ' num2str(j)];%['ReLU $5\times50$ DNN plot' num2str(j)];
    ax.Children(j).Marker = 'square';
%     ax.Children(j).MarkerSize = 8;
    ax.Children(j).LineWidth = 1.5;
%     ax.Children(j).LineStyle = '-';
    if min(ax.Children(j).XData) < xmin
        xmin = min(ax.Children(j).XData);
    end
    if max(ax.Children(j).XData) > xmax
        xmax = max(ax.Children(j).XData);
    end
    if min(ax.Children(j).YData) < ymin
        ymin = min(ax.Children(j).YData);
    end
    if max(ax.Children(j).YData) > ymax
        ymax = max(ax.Children(j).YData);
    end
end

grid off;
box on;
