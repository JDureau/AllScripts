function test = IsDefPos(M)

test = 1;
try
    if not( mean(eig(M)>0)==1 ) 
        test = 0;
    end

    n = size(M,1);
    for i = 1:n
        if 2*M(i,i)<sum(abs(M(i,:)))
            test = 0;
        end
        if 2*M(i,i)<sum(abs(M(:,i)))
            test = 0;
        end
    end
catch
    test = 0;
end