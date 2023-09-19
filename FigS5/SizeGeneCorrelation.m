clear;
clc;
close all;
CO=linspecer(3);
load('CleanChromosome.mat');
figure(1);
x=AllPlasmidSize(PlasmidIncluded==1);
y=AllPlasmidGene(PlasmidIncluded==1);
plot(x,y,'o','markersize',1,'color',CO(1,:));hold on;
set(gca,'fontsize',10);
xlabel('plasmid size','fontsize',14);
ylabel('plasmid gene number','fontsize',14);
set(gcf,'position',[100 100 220 200]);
saveas(gcf,'SizeGeneCorrelation_1.fig');
saveas(gcf,'SizeGeneCorrelation_1.pdf');


figure(2);
x=Chromosomes;
y=GeneDictionary;
plot(x(y>0),y(y>0),'o','markersize',1,'color',CO(1,:));hold on;
set(gca,'fontsize',10);
xlabel('plasmid size','fontsize',14);
ylabel('plasmid gene number','fontsize',14);
set(gcf,'position',[100 100 220 200]);
saveas(gcf,'SizeGeneCorrelation_2.fig');
saveas(gcf,'SizeGeneCorrelation_2.pdf');