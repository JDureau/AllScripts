function loss = ToOptimizeThresh(valout,valin,amplin,amplout,bootinds)

losses = [];
for indboot = 1:size(bootinds,1)
    inds = bootinds(indboot,:);
    tab = zeros(2,2);
    for i = 1:length(inds)
        if amplout(inds(i))>valout
            if amplin(inds(i))>valin
                tab(1,2) = tab(1,2)+1;
            else
                tab(1,1) = tab(1,1)+1;
            end
        else
            if amplin(inds(i))>valin
                tab(2,2) = tab(2,2)+1;
            else
                tab(2,1) = tab(2,1)+1;
            end
        end
    end
    losses(indboot) = (tab(2,2)/(tab(2,2)+tab(1,2)))^2 + (tab(1,1)/(tab(1,1)+tab(1,2)))^2;
%     losses(indboot) = (tab(2,2))^2 + (tab(1,1))^2;
%     losses(indboot) =  abs(tab(2,2)/(tab(2,2)+tab(1,2))) + abs(tab(1,1)/(tab(1,1)+tab(1,2)));
end
loss = mean(losses);
% loss = (tab(2,2))^2 + (tab(1,1))^2;


% loss = abs(tab(2,2)/(tab(2,2)+tab(1,2))) + abs(tab(1,1)/(tab(1,1)+tab(1,2)));