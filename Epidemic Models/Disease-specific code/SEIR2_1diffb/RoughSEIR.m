function Res = RoughSEIR(Data,Parameters)

NbPointsBeta = Parameters.NbPointsBeta;
LengthPointsBetaPeriod = ceil(sum(Data.NbComputingSteps)./(NbPointsBeta));
betas = ones(1,cumsum(Data.NbComputingSteps));
for i = 1:NbPointsBeta-1
    startval = Parameters.BetasValues(i);
    endval = Parameters.BetasValues(i+1);
    betas((i-1)*LengthPointsBetaPeriod+1:(i)*LengthPointsBetaPeriod) = startval:(endval-startval)/LengthPointsBetaPeriod:endval;
end


TempVariables = Parameters.InitialState;
ComputationTStep = Parameters.ComputationTStep;
TotPop = Parameters.TotalPopulation;
Variables = TempVariables;
Rebuiltf = Variables(5);
for IndTime = 2:length(Data.Instants)
    TempVariables(5) = 0;
    for IndDiscr = 1:Data.NbComputingSteps(IndTime)
        % Variables
        beta = betas(IndDiscr + sum(Data.NbComputingSteps(1:IndTime-1)));
        TempVariables(1) = TempVariables(1) + (-beta.*Variables(1).*Variables(3)/TotPop)*ComputationTStep ;% + sqrt(ComputationTStep)*Parameters.SigmaDiffusion(1)*rands(:,IndDiscr,1);
        TempVariables(2) = TempVariables(2) + ( beta.*Variables(1).*Variables(3)/TotPop-k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(2)*rands(:,IndDiscr,2);
        TempVariables(3) = TempVariables(3) + (-gamma*Variables(3) + k*Variables(2))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        TempVariables(4) = TempVariables(4) + ( gamma*Variables(3))*ComputationTStep ;%+ sqrt(ComputationTStep)*Parameters.SigmaDiffusion(3)*rands(:,IndDiscr,3);
        TempVariables(5) = TempVariables(5) + ( k*Variables(2))*ComputationTStep ;            
        Variables = TempVariables;
    end
    Rebuiltf(IndTime) = Variables(5);
end
plot(f,'g')
hold on
plot(Rebuiltf)
hold off