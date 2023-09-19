clear;
clc;
close all;
CO=linspecer(3);

load('CleanChromosome.mat');

Phis=Plasmids./Chromosomes;

data=Chromosomes(Phis<2);
ttt=Plasmids(Phis<2);

num=40;
bin=(max(data)-min(data))/num;
mea=0;
med=0;
fgh=0;
iii=fix((data-min(data))/bin);
iii(iii==num)=num-1;
iii=iii+1;
thresh=10;

for i=1:num
    xx(i)=nanmean(data(iii==i));
    mea(i)=nanmean(ttt(iii==i));
    mea_e(i)=nanstd(ttt(iii==i));
    fgh(i)=length(ttt(iii==i));
end

x=xx(fgh>thresh)';
y=log10(mea(fgh>thresh))';
plot(x,y,'.','markersize',20,'color',CO(1,:));hold on;
corr(x,y,'type','Spearman')
ft=fittype('FirstOrderHill(x, a, b)');
f=fit(x,y,ft, 'StartPoint', [2500,3]);
xf=min(x):(max(x)-min(x))/200:max(x);
yf=f(xf);
pf=predint(f,xf,0.95,'observation','off');
pfd=pf(:,1);
pfu=pf(:,2);
InBetween1=[pfd',fliplr(yf')];
InBetween2=[yf',fliplr(pfu')];
a1=fill([xf,fliplr(xf)],InBetween1,CO(1,:));hold on;
a2=fill([xf,fliplr(xf)],InBetween2,CO(1,:));hold on;
set(a1,'Facealpha',0.1);
set(a1,'Edgealpha',0);
set(a2,'Facealpha',0.1);
set(a2,'Edgealpha',0);
plot(xf,yf,'-','color',CO(1,:),'LineWidth',1);hold on;
plot(xf,pfd,'--','color',CO(1,:),'LineWidth',1);hold on;
plot(xf,pfu,'--','color',CO(1,:),'LineWidth',1);hold on;



xlabel('host chromosome size','fontsize',14);
ylabel('plasmid contents','fontsize',14);
legend off;
set(gcf,'position',[100 100 210 200]);
saveas(gcf,'PlasmidContent_Bin40.fig');
saveas(gcf,'PlasmidContent_Bin40.pdf');




