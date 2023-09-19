clear;
clc;
close all;
num_meta=80;% number of subtances

NUMGENEs=[100 500 1000]-num_meta/4;
REP=20;
Plas=0:50:500;
DIS=0*ones(length(NUMGENEs),REP,length(Plas));
BURD=0*ones(length(NUMGENEs),REP,length(Plas));

for sdf=1:length(NUMGENEs)
    
    NUMGENE=NUMGENEs(sdf);
    genome=NUMGENE+num_meta/4;% number of genes
    for rty=1:REP
        (sdf-1)*REP+rty
        S_mother=GenomeMatrix(num_meta,genome);
        num_fluxs=genome+Plas;
        coef0=rand(1,max(num_fluxs));
        GR=0*num_fluxs; 
        FFG=0*ones(length(num_fluxs),genome);
        for i=1:length(num_fluxs)
            num_flux=num_fluxs(i);
            coef=[coef0(1:genome) 5*coef0(genome+1:num_flux)];

            % clear x,S,con,prob,con;
            prob=optimproblem('ObjectiveSense','max');
        
            x = optimvar('x',1,num_flux,'LowerBound',0,'UpperBound',1);
        
            S=[S_mother PlasmidMatrix(num_meta,num_flux-genome)];
            wer=sum(times(S,repelem(x,num_meta,1)),2);
            con=wer(1:num_meta/2)==0*ones(num_meta/2,1);
            prob.Constraints.con=con;
    
            NUMGENE=genome-num_meta/4;    
            prob.Objective=sum(coef(3/4*NUMGENE+1+num_meta/4:num_flux).*x(3/4*NUMGENE+1+num_meta/4:num_flux));
            % prob.Objective=sum(coef(3/4*NUMGENE+1+num_meta/4:num_flux).*x(3/4*NUMGENE+1+num_meta/4:num_flux));
            x0.x=rand(1,num_flux);
            sol=solve(prob,x0);
            tt=sol.x;
            FFG(i,:)=tt(1,1:genome);
            if ~isempty(tt)
                GR(i)=sum(coef(3/4*NUMGENE+num_meta/4+1:genome).*tt(3/4*NUMGENE+num_meta/4+1:genome));
            end
            DIS(sdf,rty,i)=BrayCurtis(FFG(1,:),FFG(i,:));
            BURD(sdf,rty,i)=GR(i)/GR(1);
        end
        
    end
end
CCC=linspecer(length(NUMGENEs));

for i=1:length(NUMGENEs)
    mm=mean(BURD(i,:,:),2);
    rr=std(BURD(i,:,:),0,2);
    errorbar(Plas,reshape(mm,1,length(Plas)),reshape(rr,1,length(Plas)),'o-','linewidth',2,'markersize',8,'color',CCC(i,:),'CapSize',0);hold on;
end

xlim([1,max(Plas)]);
set(gca,'fontsize',16);
xlabel('plasmid size','fontsize',20);
ylabel('growth rate','fontsize',20);
set(gcf,'position',[100 100 300 300]);
saveas(gcf,'BurdenPlasmidSize.fig');
exportgraphics(gcf,'BurdenPlasmidSize.pdf');
