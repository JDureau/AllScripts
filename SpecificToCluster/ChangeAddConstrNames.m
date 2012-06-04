
for i = 1:100
    try
        load(['/VIH_PlayingWithSigm_ResultsAddConstr__1_' num2str(i) '.mat'])
        save(['/VIH_PlayingWithSigm_ResultsSigmoid__1_' num2str(i) '.mat'],'Res')
    
    end
end