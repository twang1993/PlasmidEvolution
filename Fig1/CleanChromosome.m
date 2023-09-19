clear;
clc;
load('DataSets.mat');
AllPlasmidSize=T1{:,8};
AllPlasmidGC=T1{:,9}/100.*AllPlasmidSize;
AllPlasmidProtein=str2double(T1{:,10});
AllPlasmidGene=str2double(T1{:,14});
AllPlasmidHostName=T1{:,1};
AllPlasmidRefSeq=T1{:,6};
AllPlasmidINSDC=T1{:,7};

OrganismNameDictionary1=T2{:,1};
RepliconDictionary1=T2{:,10};
SizeDictionary1=T2{:,8}*10^3;
GCDictionary1=T2{:,9}/100.*SizeDictionary1;
GeneDictionary1=T2{:,17};
GroupDictionary1=T2{:,2};
KingdomDictionary1=T2{:,2};
PhylumDictionary1=T2{:,2};

for i=1:length(KingdomDictionary1)
    i
    temp=GroupDictionary1{i};
    c=strfind(temp,';');
    KingdomDictionary1{i}=temp(1:c(1)-1);
    PhylumDictionary1{i}=temp(c(1)+1:c(2)-1);
    SubPhylumDictionary1{i}=temp(c(2)+1:end);
end

IsComplete=contains(RepliconDictionary1,'chromosome').*contains(RepliconDictionary1,'plasmid');
% MultiChromosome=contains(RepliconDictionary1,'chromosome 2');

OrganismNameDictionary1=OrganismNameDictionary1(IsComplete>0);
RepliconDictionary1=RepliconDictionary1(IsComplete>0);
KingdomDictionary1=KingdomDictionary1(IsComplete>0);
PhylumDictionary1=PhylumDictionary1(IsComplete>0);
SubPhylumDictionary1=SubPhylumDictionary1(IsComplete>0);
GroupDictionary1=GroupDictionary1(IsComplete>0);
SizeDictionary1=SizeDictionary1(IsComplete>0);
GeneDictionary1=GeneDictionary1(IsComplete>0);
GCDictionary1=GCDictionary1(IsComplete>0);

GenomePlasmidNumber1=0*SizeDictionary1;

for i=1:length(GenomePlasmidNumber1)
    i
    temp=string(RepliconDictionary1(i));
    c=strfind(temp,'; plasmid');
    if isempty(c)==0
        GenomePlasmidNumber1(i)=length(c);
    end
end

GenomePlasmidNumber=GenomePlasmidNumber1;
SizeDictionary=SizeDictionary1;
GeneDictionary=GeneDictionary1;
GCDictionary=GCDictionary1;
PlasmidIncluded=0*ones(length(AllPlasmidSize),1);

pin=1;
PPP=1:length(AllPlasmidSize);
PCmap=0;

for i=1:length(AllPlasmidSize)
    i
    temp1=string(AllPlasmidRefSeq(i));
    temp2=string(AllPlasmidINSDC(i));
    c1=0*ones(length(SizeDictionary),1);
    c2=c1;
    if length(char(temp1))>3
        c1=contains(RepliconDictionary1,temp1);
    end
    if length(char(temp2))>3
        c2=contains(RepliconDictionary1,temp2);
    end
    
    if sum(c1+c2)>0
        PlasmidIncluded(i)=1;
    end
        GenomePlasmidNumber=GenomePlasmidNumber-(c1+c2>0);
        if ~isnan(AllPlasmidSize(i))
            SizeDictionary=SizeDictionary-(c1+c2>0)*AllPlasmidSize(i);
        end
        if ~isnan(AllPlasmidGene(i))
            GeneDictionary=GeneDictionary-(c1+c2>0)*AllPlasmidGene(i);
        end
        if ~isnan(AllPlasmidGC(i))
            GCDictionary=GCDictionary-(c1+c2>0)*AllPlasmidGC(i);
        end
      
     temp=PPP(c1+c2>0);
     for cvb=1:length(temp)
        PCmap(1,pin)=temp(cvb);
        PCmap(2,pin)=i;
        pin=pin+1;
     end
end

T3=readtable('ChromosomeManual_JKH.xlsx','VariableNamingRule','preserve');
ManualIndex=T3{:,1};
SizeDictionary(ManualIndex)=T3{:,4}/10^3;
GenomePlasmidNumber(ManualIndex)=0;

Index=GenomePlasmidNumber==0;

Chromosomes=SizeDictionary(Index);
Genomes=SizeDictionary1(Index);
Plasmids=Genomes-Chromosomes;
GCs=GCDictionary(Index)./Chromosomes;
KingdomDictionary=KingdomDictionary1(Index);
PhylumDictionary=PhylumDictionary1(Index);
SubPhylumDictionary=SubPhylumDictionary1(Index);
PlasmidNum=GenomePlasmidNumber1(Index);

save('CleanChromosome.mat');

