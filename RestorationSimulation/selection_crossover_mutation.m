function newchrom = selection_crossover_mutation(chrom,fitnv,Npg,selp,crossp,mutp)
[Nc,~] = size(chrom);
newchrom = chrom*0;
fprob = cumsum(fitnv/sum(fitnv)); % RandStream.setDefaultStream(RandStream('mt19937ar','Seed',sum(1000*clock)));
%selection according to the fitness value
for c = 1 : ceil(Nc*selp)
    [~,fi] = max(fitnv);
    newchrom(c,:)=chrom(fi,:);
    fitnv(fi)=0;
end
c = ceil(Nc*selp)+1;
while c <= Nc
    temp = find(fprob>=rand());
    newchrom(c,:) = chrom(temp(1),:);
    c=c+1;
end
% crossover
c = ceil(Nc*selp)+1;
while c+1 <= Nc
    if Npg(1) > 0
        chrom1 = newchrom(c,1:Npg(1));
        chrom2 = newchrom(c+1,1:Npg(1));
        if rand() < crossp
            pst = ceil(rand()*(Npg(1)-1));
            if pst < (Npg(1)-1)/2
                remain1 = chrom1(pst+1:Npg(1));remain2=chrom2(pst+1:Npg(1));
                new1 = getdiff(chrom2,remain1);
                new2 = getdiff(chrom1,remain2);
                newchrom1 = [new1 remain1];
                newchrom2 = [new2 remain2];
            else
                remain1 = chrom1(1:pst);
                remain2 = chrom2(1:pst);
                new1 = getdiff(chrom2,remain1);
                new2 = getdiff(chrom1,remain2);
                newchrom1 = [remain1 new1];
                newchrom2 = [remain2 new2];
            end
        else
            newchrom1 = chrom1;
            newchrom2 = chrom2;
        end
    else
        newchrom1 = [];
        newchrom2 = [];
    end

    if Npg(2) > 0
        chrom3 = newchrom(c,Npg(1)+1:Npg(1)+Npg(2));
        chrom4 = newchrom(c+1,Npg(1)+1:Npg(1)+Npg(2));
        if rand() < crossp
            pst = ceil(rand()*(Npg(2)-1));
            if pst < (Npg(2)-1)/2
                remain1 = chrom3(pst+1:Npg(2));
                remain2 = chrom4(pst+1:Npg(2));
                new1 = getdiff(chrom4,remain1);new2=getdiff(chrom3,remain2);
                newchrom3 = [new1 remain1];
                newchrom4 = [new2 remain2];
            else
                remain1 = chrom3(1:pst);
                remain2 = chrom4(1:pst);
                new1 = getdiff(chrom4,remain1);
                new2 = getdiff(chrom3,remain2);
                newchrom3 = [remain1 new1];
                newchrom4 = [remain2 new2];
            end
        else
            newchrom3 = chrom3;
            newchrom4 = chrom4;
        end
    else
        newchrom3 = [];
        newchrom4 = [];
    end
    newchrom(c,:) = [newchrom1 newchrom3];
    newchrom(c+1,:) = [newchrom2 newchrom4];
    c = c+2;
end
% mutation
for c = 1 : Nc
    if c > ceil(Nc*selp)
        if rand() < mutp && Npg(1) >= 2
            newseq = randperm(Npg(1));
            temp = newchrom(c,newseq(1));
            newchrom(c,newseq(1)) = newchrom(c,newseq(2));
            newchrom(c,newseq(2)) = temp;
        end
        if rand() < mutp && Npg(2) >= 2
            newseq = randperm(Npg(2));
            temp = newchrom(c,Npg(1)+newseq(1));
            newchrom(c,Npg(1)+newseq(1)) = newchrom(c,Npg(1)+newseq(2));
            newchrom(c,Npg(1)+newseq(2)) = temp;
        end
    end
end
end

function vectorA = getdiff(vectorA,vectorB)
[~,op] = intersect(vectorA,vectorB);
vectorA(op) = [];
end
