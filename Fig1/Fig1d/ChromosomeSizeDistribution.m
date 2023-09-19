clear;
clc;
close all;

load('CleanChromosome.mat');

histogram(log10(Chromosomes),150,'Normalization','probability','EdgeAlpha',0);
xlabel('chromosome size','fontsize',14);
ylabel('frequency','fontsize',14);
xlim([2.5 4.2]);
set(gcf,'position',[10 10 200 200]);
saveas(gcf,'ChromosomeSizeDistribution.fig');
exportgraphics(gcf,'ChromosomeSizeDistribution.pdf');