clear;
clc;
close all;

explode1 = [0 1];
label1=[{'no plasmids'},{'at least one plasmid'}];
pie([36593 15713],explode1);
ax=gca,
for i=1:2
        ax.Children(2*i).EdgeAlpha = 0;
end
delete(findobj(ax,'Type','text'));
colormap lines;
legend(label1,'location','eastoutside');
legend box off;
saveas(gcf,'WithPlasmidPie.fig');
exportgraphics(gcf,'WithPlasmidPie.pdf');