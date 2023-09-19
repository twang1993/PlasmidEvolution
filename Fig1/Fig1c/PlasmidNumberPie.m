clear;
clc;
close all;
load('CleanChromosome.mat');

PlasmidNum(PlasmidNum>5)=6;
PlsN=unique(PlasmidNum);
PlsNIndex=0*ones(size(PlasmidNum,1),1);
for i=1:size(PlsN,1)
    c=find(ismember(PlasmidNum,PlsN(i)));
    PlsNIndex(c)=i;
    PN(i)=length(c);
end

label1=num2str(PlsN);
pie(PN);
ax=gca,
for i=1:length(PN)
        ax.Children(2*i).EdgeAlpha = 0;
end
ax.Colormap=linspecer(length(PN)+1);
delete(findobj(ax,'Type','text'));
legend(label1,'location','eastoutside');
legend box off;
saveas(gcf,'PlasmidNumberPie.fig');
exportgraphics(gcf,'PlasmidNumberPie.pdf');