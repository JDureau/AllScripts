function loss = ToOptimizeThreshInfl(propout,propin,propsin,propsout,amplin,valampl)

tab = zeros(2,2);
for i = 1:length(amplin)
    if amplin(i)>valampl
        if propsout(i)>propout
            if propsin(i)>propin
                tab(1,2) = tab(1,2)+1;
            else
                tab(1,1) = tab(1,1)+1;
            end
        else
            if propsin(i)>propin
                tab(2,2) = tab(2,2)+1;
            else
                tab(2,1) = tab(2,1)+1;
            end
        end
    end
end

loss = (tab(2,2)/(tab(2,2)+tab(1,2)))^2 + (tab(1,1)/(tab(1,1)+tab(1,2)))^2;