function [] = RunHIVSimulations(BigInd,SmallInd)




load([pwd '/AllData/Avahan/VIH_PlayingWithSigm_ResGens' num2str(BigInd) '_Biased.mat'])

inds = 10*(SmallInd-1)+1:10*(SmallInd);
Ress = {};
for i = inds(1):inds(end)
 
%       if BigInd == 1
% %	BiasedParameters = UpdateParsNoTransfToTransf(ResGens1{i}.BiasedParameters);
% 	Res = HIVapplyInference(ResGens1{i}.Data,ResGens1{i}.BiasedParameters);
%       else
% %	BiasedParameters = UpdateParsNoTransfToTransf(ResGens2{i}.BiasedParameters);
% 	Res = HIVapplyInference(ResGens2{i}.Data,ResGens2{i}.BiasedParameters);
%       end
%         Res = HIVapplyInference(ResGen.Data,ParametersBelgaum);
% 
    Res = 10;
    Ress{end+1} = Res;
    %         
    %         
    save([pwd '/AllData/Avahan/VIH_PlayingWithSigm_RightPriors' num2str(BigInd) '_' num2str(SmallInd) '_Biased.mat'],'Ress')
end