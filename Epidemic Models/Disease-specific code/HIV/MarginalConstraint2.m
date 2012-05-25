function loss = MarginalConstraint2(valout,valin,amplin,amplout,objective,bootinds)



gooddetprops = [];
for indboot = 1:size(bootinds,1)
    tab = zeros(2,2);
    for i = 1:size(bootinds,2)
        ind = bootinds(indboot,i);
        if amplout(ind)>valout
            if amplin(ind)>valin
                tab(1,2) = tab(1,2)+1;
            else
                tab(1,1) = tab(1,1)+1;
            end
        else
            if amplin(ind)>valin
                tab(2,2) = tab(2,2)+1;
            else
                tab(2,1) = tab(2,1)+1;
            end
        end
    end
    gooddetprops(indboot) = (tab(1,2)/(tab(1,1)+tab(1,2)));

end
% tab
% mean(gooddetprops)
if mean(gooddetprops)<objective
    loss = 10*abs(mean(gooddetprops)-objective);
%     loss
else
    loss = abs(mean(gooddetprops)-objective);
end