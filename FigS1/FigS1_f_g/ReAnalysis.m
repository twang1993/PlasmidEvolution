clear;
clc;
close all;
load('CleanChromosome.mat');

Phylums=unique(PhylumDictionary);
PhylumIndex=0*ones(size(PhylumDictionary,1),1);
for i=1:size(Phylums,1)
    c=find(ismember(PhylumDictionary,Phylums(i)));
    PhylumIndex(c)=i;
    PN(i)=length(c);
end


Phis=Plasmids./Chromosomes;


CO=linspecer(length(Phylums));

figure(1);
x=0;
y=0;
for i=1:length(Phylums)
    y(i)=mean(Plasmids(PhylumIndex==i&Phis<2));
    y_err(i)=std(Plasmids(PhylumIndex==i&Phis<2));
    x(i)=mean(Chromosomes(PhylumIndex==i));
    x_err(i)=std(Chromosomes(PhylumIndex==i));
    errorbar(x(i),y(i),y_err(i),'o','markersize',2*max(log(PN(i)),0.1),'CapSize',0,'linewidth',1.5,'color',CO(i,:));hold on;

end
set(gca,'YScale','log');
xlim([1000 8000]);
ylim([10^1 10^3.2]);
set(gca,'FontSize',10);
xlabel('mean chromosome size','FontSize',14);
ylabel('total plasmid size','FontSize',14);
set(gcf,'position',[100 100 300 300]);
saveas(gcf,'ReAnalysis_1.fig');
saveas(gcf,'ReAnalysis_1.pdf');

figure(2);
ttt=Plasmids./PlasmidNum;
x=0;
y=0;
for i=1:length(Phylums)
    y(i)=mean(ttt(PhylumIndex==i&Phis<2));
    y_err(i)=std(ttt(PhylumIndex==i&Phis<2));
    x(i)=mean(Chromosomes(PhylumIndex==i));
    x_err(i)=std(Chromosomes(PhylumIndex==i));
    errorbar(x(i),y(i),y_err(i),'o','markersize',2*max(log(PN(i)),0.1),'CapSize',0,'linewidth',1.5,'color',CO(i,:));hold on;

end

xlim([1000 8000]);
ylim([10^0.6 10^2.8]);
set(gca,'YScale','log');
set(gca,'FontSize',10);
xlabel('mean chromosome size','FontSize',14);
ylabel('mean size per plasmid','FontSize',14);
set(gcf,'position',[100 100 300 300]);
saveas(gcf,'ReAnalysis_2.fig');
saveas(gcf,'ReAnalysis_2.pdf');

