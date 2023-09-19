clear;
clc;
close all;
global mu kappa D;
mu=0.5;
kappas=[0.001 0.005 0.01 0.05];
D=0.1;
initial=[0,1];
timespan=0:0.1:1000;
CO=linspecer(length(kappas));
for i=1:length(kappas)
    kappa=kappas(i);
    [t,y]=ode45(@SegregationError,timespan,initial);
    plot(t,y(:,2)./(y(:,1)+y(:,2)),'.','Color',CO(i,:),'markersize',6);hold on;
end
set(gca,'Fontsize',10);
xlabel('time','fontsize',14);
ylabel('plasmid abundance','fontsize',14);
set(gcf,'position',[100 100 250 220]);
saveas(gcf,'PlasmidLoss.fig');
saveas(gcf,'PlasmidLoss.pdf');
exportgraphics(gcf,'PlasmidLoss.jpg','Resolution',1000);


function dydt=SegregationError(t,y)
global mu kappa D;
dydt=[mu*y(1)*(1-y(1)-y(2))+kappa*y(2)-D*y(1);
    mu*y(2)*(1-y(1)-y(2))-kappa*y(2)-D*y(2)];
end