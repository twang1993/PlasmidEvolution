function S=PlasmidMatrix(NumMet,NumGene)
    S=0*ones(NumMet,NumGene);
    range=3;
    Num1=NumMet/4;

    t=S(1:2*Num1,:);
    S(1:2*Num1,:)=spar(t,0.3,-1,3);

    t=S(Num1*3+1:4*Num1,:);
    S(Num1*3+1:4*Num1,:)=spar(t,0.3,1,3);

end