clear;
clc;
close all;

load('CleanChromosome.mat');

histogram(log10(AllPlasmidSize(PlasmidIncluded>0)),100,'Normalization','probability','EdgeAlpha',0);
xlabel('plasmid size','fontsize',14);
ylabel('frequency','fontsize',14);
xlim([0 4]);
set(gcf,'position',[10 10 200 200]);
saveas(gcf,'PlasmidSizeDistribution.fig');
exportgraphics(gcf,'PlasmidSizeDistribution.pdf');