function y=spar(x,f,pos,range)
    y=0*x;
    for i=1:size(x,2)
        temp=0*ones(size(x,1),1);
        ii=randperm(size(x,1));
        nn=max(1,fix(rand*size(x,1)*f));
        if pos==1
            for j=1:nn
                y(ii(j),i)=rndin(1:range);
            end
        end
        if pos==-1
            for j=1:nn
                y(ii(j),i)=-rndin(1:range);
            end
        end
        if pos==0
            for j=1:nn
                y(ii(j),i)=rndin(1:range)*((rand>0.5)*2-1);
            end
        end
    end
end