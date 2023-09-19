clear;
clc;
close all;
load('CleanChromosome.mat');

xdata=0*ones(size(PCmap,2),1);
ydata=0*ones(size(PCmap,2),1);
zdata=0*ones(size(PCmap,2),1);
for i=1:size(PCmap,2)
    if Index(PCmap(1,i))==1
        xdata(i)=SizeDictionary(PCmap(1,i));
        zdata(i)=SizeDictionary1(PCmap(1,i))-SizeDictionary(PCmap(1,i));
        ydata(i)=AllPlasmidSize(PCmap(2,i));
    end
end

data=xdata(zdata./xdata<2);
ttt=ydata(zdata./xdata<2);
num=10;
bin=(max(data)-min(data))/num;

iii=fix((data-min(data))/bin);
iii(iii==num)=num-1;
iii=iii+1;
thresh=50;
xx=0;
mea=0;
mea_e=0;
fgh=0;
frc10=0;
frc20=0;
frc30=0;
frc40=0;
frc50=0;

for i=1:num
    xx(i)=nanmean(data(iii==i));
    mea(i)=nanmean(ttt(iii==i));
    mea_e(i)=nanstd(ttt(iii==i));
    fgh(i)=length(ttt(iii==i));
    frc8(i)=sum(ttt(iii==i)<8);
    frc10(i)=sum(ttt(iii==i)<10);
    frc20(i)=sum(ttt(iii==i)<20);
    frc50(i)=sum(ttt(iii==i)<50);
end


figure(1);
for i=1:num
    subplot(num,1,i);
    if fgh(i)>thresh
    histogram(log10(ttt(iii==i)),[0:0.1:3],'Normalization','probability','EdgeAlpha',0);hold on;
    end
    xticks([0:0.5:3]);
    xticklabels([]);
    yticks([]);
    ylim([0 0.2]);
    yticks([0 0.1 0.2]);
    yticklabels([]);
end
set(gcf,'position',[100 100 210 400]);
saveas(gcf,'SmallPlasmidFraction_1.fig');
saveas(gcf,'SmallPlasmidFraction_1.pdf');

figure(2);
CO=linspecer(4);
plot(xx(fgh>thresh),frc8(fgh>thresh)./fgh(fgh>thresh),'o-','color',CO(1,:),'linewidth',1,'markersize',8);hold on;
plot(xx(fgh>thresh),frc10(fgh>thresh)./fgh(fgh>thresh),'o-','color',CO(2,:),'linewidth',1,'markersize',8);hold on;
plot(xx(fgh>thresh),frc20(fgh>thresh)./fgh(fgh>thresh),'o-','color',CO(3,:),'linewidth',1,'markersize',8);hold on;
plot(xx(fgh>thresh),frc50(fgh>thresh)./fgh(fgh>thresh),'o-','color',CO(4,:),'linewidth',1,'markersize',8);hold on;
set(gca,'fontsize',12);
box off;
legend('8','10','20','50');
legend box off;
xlabel('chromosome size','fontsize',14);
ylabel('fraction of small plasmids','fontsize',14);
set(gcf,'position',[100 100 250 400]);
saveas(gcf,'SmallPlasmidFraction_2.fig');
saveas(gcf,'SmallPlasmidFraction_2.pdf');