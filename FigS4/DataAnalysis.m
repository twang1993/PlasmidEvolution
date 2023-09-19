clear;
clc;
close all;
T=readtable('exoplanet_eu_catalog.csv');

StarName=T{:,69};
StarMass=T{:,83};
% StarMass=T{:,86};
PlanetMass=T{:,3};

StarUniqueName=unique(StarName);
NumPlanetPerStar=0*ones(length(StarUniqueName),1);
StarMassUnique=0*ones(length(StarUniqueName),1);
PlanetMassPerStar=0*ones(length(StarUniqueName),1);
MassPerPlanet=0*ones(length(StarUniqueName),1);

for i=1:size(StarUniqueName)
    i
    temp=strcmp(StarName,StarUniqueName{i});
    NumPlanetPerStar(i)=sum(temp);
    StarMassUnique(i)=mean(StarMass(temp));
    PlanetMassPerStar(i)=sum(PlanetMass(temp));
    MassPerPlanet(i)=PlanetMassPerStar(i)/NumPlanetPerStar(i);
end
figure(1);
plot(StarMassUnique(2:end),NumPlanetPerStar(2:end),'ko','markersize',10);
set(gca,'XScale','log');
set(gca,'fontsize',12);
axis([10^(-2) 10^1.5 0 8]);
xlabel('Star Mass','fontsize',14);
ylabel('number of oribiting planets','fontsize',14);
set(gcf,'position',[100 100 800 200]);

saveas(gcf,'DataAnalysis.fig');
saveas(gcf,'DataAnalysis.pdf');

%%%%%%%%%%%%%%%%%
data=StarMassUnique(2:end);
ttt=PlanetMassPerStar(2:end);
plot(data,ttt,'k.','markersize',5);hold on;
corr(data(data.*ttt>0),ttt(data.*ttt>0),'type','Spearman')
set(gca,'XScale','log');
ylim([0 20]);
xlim([10^(-1.2) 10^0.5]);
set(gca,'fontsize',8);
xlabel('Star Mass','fontsize',14);
ylabel('total planet mass','fontsize',14);
set(gcf,'position',[100 100 350 200]);

saveas(gcf,'DataAnalysis_2.fig');
saveas(gcf,'DataAnalysis_2.pdf');

figure(3);
data=StarMassUnique(2:end);
ttt=MassPerPlanet(2:end);
plot(data,ttt,'k.','markersize',5);
corr(data(data.*ttt>0),ttt(data.*ttt>0),'type','Spearman')
set(gca,'fontsize',8);
set(gca,'XScale','log');
ylim([0 20]);
xlim([10^(-1.2) 10^0.5]);
xlabel('Star Mass','fontsize',14);
ylabel('mean planet mass','fontsize',14);
set(gcf,'position',[100 100 350 200]);

saveas(gcf,'DataAnalysis_3.fig');
saveas(gcf,'DataAnalysis_3.pdf');
