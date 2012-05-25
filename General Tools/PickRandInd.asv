function ind = PickRandInd(Weigths)

if sum(Weigths)>0
    Weigths = Weigths/sum(Weigths);
else
    Weigths = 1/length(Weigths)*ones(size(Weigths));
end
    
rd = rand(1,1);
test = 0;
ind = 1;
while not(test)
    if rd < sum(Weigths(1:ind))
        test = 1;
    else
        ind = ind+1;
    end
end
