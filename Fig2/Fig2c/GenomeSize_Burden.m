clear;
clc;
close all;
num_meta=80;% number of subtances

NUMGENEs=[100:200:10000]-num_meta/4;
REP=20;
ratios=[1 1.2];
DIS=0*ones(length(NUMGENEs),REP,length(ratios));
BURD=0*ones(length(NUMGENEs),REP,length(ratios));

for sdf=1:length(NUMGENEs)
    sdf
    NUMGENE=NUMGENEs(sdf);
    genome=NUMGENE+num_meta/4;% number of genes
    for rty=1:REP
        S_mother=GenomeMatrix(num_meta,genome);
        num_fluxs=fix(genome*ratios);
        coef0=rand(1,max(num_fluxs));
        GR=0*num_fluxs; 
        FFG=0*ones(length(num_fluxs),genome);
        for i=1:length(num_fluxs)
            i
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
figure(1);

mm=mean(BURD(:,:,2),2);
rr=std(BURD(:,:,2),0,2);
errorbar(NUMGENEs(1:2:end),mm(1:2:end),rr(1:2:end),'o-','linewidth',2,'markersize',6,'color',CCC(i,:),'CapSize',0);hold on;

xlim([1,max(NUMGENEs)]);
set(gca,'fontsize',16);
xlabel('genome size','fontsize',20);
ylabel('growth rate','fontsize',20);
set(gcf,'position',[100 100 300 300]);
saveas(gcf,'GenomeSize_Burden.fig');
exportgraphics(gcf,'GenomeSize_Burden.pdf');
