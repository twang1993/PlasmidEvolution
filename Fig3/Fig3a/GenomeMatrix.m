function S=GenomeMatrix(NumMet,NumGene)
    NumGeneIn=NumGene-NumMet/4;
    S=0*ones(NumMet,NumGeneIn+NumMet/4);
    range=3;
    Num1=NumMet/4;

    t=S(1:Num1,1+NumMet/4:NumGeneIn/4+NumMet/4);
    S(1:Num1,1+NumMet/4:NumGeneIn/4+NumMet/4)=spar(t,0.3,-1,3);

    t=S(Num1+1:2*Num1,1+NumMet/4:NumGeneIn/4+NumMet/4);
    S(Num1+1:2*Num1,1+NumMet/4:NumGeneIn/4+NumMet/4)=spar(t,0.3,1,3);

    t=S(Num1+1:2*Num1,NumGeneIn/4+1+NumMet/4:NumGeneIn*3/4+NumMet/4);
    S(Num1+1:2*Num1,NumGeneIn/4+NumMet/4+1:NumGeneIn*3/4+NumMet/4)=spar(t,0.3,0,3);

    t=S(Num1+1:2*Num1,NumGeneIn*3/4+1+NumMet/4:NumGeneIn+NumMet/4);
    S(Num1+1:2*Num1,NumGeneIn*3/4+1+NumMet/4:NumGeneIn+NumMet/4)=spar(t,0.3,-1,3);

    t=S(2*Num1+1:3*Num1,NumGeneIn*3/4+1+NumMet/4:NumGeneIn+NumMet/4);
    S(2*Num1+1:3*Num1,NumGeneIn*3/4+NumMet/4+1:NumGeneIn+NumMet/4)=spar(t,0.3,1,3);

    for i=1:NumMet/4
        S(i,i)=rndin([1:3]);
    end
end

