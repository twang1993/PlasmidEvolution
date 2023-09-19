clear;
clc;
close all;
num_meta=80;

NUMGENEs=[100:500:15000]-num_meta/4;
qws=2.^[0:1:10];
Plas=fix([0 qws]);
Rep=10;
PlasIndex=0*ones(length(NUMGENEs),length(Plas),Rep);

for ert=1:Rep
    for sdf=1:length(NUMGENEs)
        NUMGENE=NUMGENEs(sdf);
        genome=NUMGENE+num_meta/4;
    
        S_mother=GenomeMatrix(num_meta,genome);
        GR=0*ones(length(Plas),1);
        for i=1:length(Plas)
            i+(sdf-1)*length(Plas)
            if i==1
                num_flux=genome+Plas(i);
            else
                num_flux=genome+sum(Plas(1:i-1).*PlasIndex(sdf,1:i-1,ert))+Plas(i);
            end
            coef0=rand(1,num_flux);
            coef=[coef0(1:genome) 5*coef0(genome+1:num_flux)];
            prob=optimproblem('ObjectiveSense','max');
            x = optimvar('x',1,num_flux,'LowerBound',0,'UpperBound',1);
            S=[S_mother PlasmidMatrix(num_meta,num_flux-genome)];
            wer=sum(times(S,repelem(x,num_meta,1)),2);
            con=wer(1:num_meta/2)==0*ones(num_meta/2,1);
            prob.Constraints.con=con;
    
            NUMGENE=genome-num_meta/4;    
            prob.Objective=sum(coef(3/4*NUMGENE+1+num_meta/4:num_flux).*x(3/4*NUMGENE+1+num_meta/4:num_flux));
            x0.x=rand(1,num_flux);
            sol=solve(prob,x0);
            tt=sol.x;
            if ~isempty(tt)
                ChroIndex=3/4*NUMGENE+num_meta/4+1:genome;
                GR(i)=sum(coef(ChroIndex).*tt(ChroIndex));
                if i==1
                    PGR=GR(i);
                else
                    PGR=sum(coef(num_flux-Plas(i)+1:num_flux).*tt(num_flux-Plas(i)+1:num_flux));
                end
            end
            Burd=GR(i)/GR(1);
            if (Burd>0.8)&&(PGR/GR(i)>0.05)
                PlasIndex(sdf,i,ert)=1;
            end

            if i>=2&PlasIndex(sdf,i,ert)==1
                for j=2:i
                    if PlasIndex(sdf,j,ert)==1
                        break;
                    end
                end
                sm=0;
                for k=j:i
                    if PlasIndex(sdf,k,ert)==1
                        PGR=sum(coef(genome+sm+1:genome+sm+Plas(k)).*tt(genome+sm+1:genome+sm+Plas(k)));
                        if PGR/GR(i)<0.05
                            PlasIndex(sdf,k,ert)=0;
                        end
                        sm=sm+Plas(k);
                    end
                end

            end


            % if i>=3
            %     if (PlasIndex(sdf,i,ert)==0)&&(PlasIndex(sdf,i-1,ert)==1)
            %         break;
            %     end
            % end

        end
    end
end

PlasContent=0*ones(length(NUMGENEs),Rep);
for i=1:length(NUMGENEs)
    for j=1:Rep
        PlasContent(i,j)=sum(PlasIndex(i,:,j).*Plas);
    end
end

PlasNumber=0*ones(length(NUMGENEs),Rep);
for i=1:length(NUMGENEs)
    for j=1:Rep
        PlasNumber(i,j)=nnz(PlasIndex(i,:,j))-1;
    end
end

CO=linspecer(3);
figure(1);
errorbar(NUMGENEs+num_meta/4,mean(PlasContent,2),std(PlasContent,0,2),'.-','color',CO(1,:),'CapSize',0,'MarkerSize',20,'LineWidth',1);hold on;
set(gca,'YScale','log');
xlabel('host chromosome size','fontsize',14);
ylabel('total plasmid size','fontsize',14);
set(gcf,'position',[100 100 210 200]);
saveas(gcf,'CombinedAll_1.fig');
saveas(gcf,'CombinedAll_1.pdf');

figure(2);
errorbar(NUMGENEs+num_meta/4,mean(PlasNumber,2),std(PlasNumber,0,2),'.-','color',CO(2,:),'CapSize',0,'MarkerSize',20,'LineWidth',1);hold on;
xlabel('host chromosome size','fontsize',14);
ylabel('mean plasmid number','fontsize',14);
set(gcf,'position',[100 100 210 200]);
saveas(gcf,'CombinedAll_2.fig');
saveas(gcf,'CombinedAll_2.pdf');

figure(3);
errorbar(NUMGENEs+num_meta/4,nanmean(PlasContent./PlasNumber,2),nanstd(PlasContent./PlasNumber,0,2),'.-','color',CO(3,:),'CapSize',0,'MarkerSize',20,'LineWidth',1);hold on;
set(gca,'YScale','log');
xlabel('host chromosome size','fontsize',14);
ylabel('mean plasmid size','fontsize',14);
set(gcf,'position',[100 100 210 200]);
saveas(gcf,'CombinedAll_3.fig');
saveas(gcf,'CombinedAll_3.pdf');


