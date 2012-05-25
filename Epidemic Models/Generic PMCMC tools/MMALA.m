function res = MMALA(Data, Parameters, epsil)

N = size(Data.Observations,2);

% Initialise TempPar
ResTemp = EstimationSMCsmoothAndGrad(Data,Parameters);
TempPar.LogLik = ResTemp.LogLik;
GradTemp = ResTemp.Grad;
if Parameters.DiffusionType == 'Add'
    if Parameters.EstimatedVariables.SigmaRW 
        TempPar.SigmaRW = Parameters.SigmaRW;
        GTemp = (2*N)/TempPar.SigmaRW^2;
        TempPar.Mu = TempPar.SigmaRW + Parameters.Epsil^2*GTemp^-1*GradTemp/2 - epsil^2*SigmaRWTemp/2;
        TempPar.Sigma = epsil*sqrt(GTemp^-1);
    end
end

% Calibrate the method
[Parameters, TempPar] = CalibrateMethod( Data, Parameters, TempPar);

% Burn-In



Thetas = [];
LogLiks = [];

AccRate = 1;
Condition = 0;
while not(Condition)
    Accepted = [];
    for IndIt = 1:20
        if Parameters.DiffusionType == 'Add'
            if Parameters.EstimatedVariables.SigmaRW 
                SigmaRWStar = MuTemp + SigmaTemp*randn(1,1);
                ParametersStar = ParametersTemp;
                ParametersStar.SigmaRW = SigmaRWStar;
                GStar = (2*N)/SigmaRWStar^2;
            end
        end
        ResStar = EstimationSMCsmoothAndGrad(Data,ParametersStar);
        LogLikStar = ResStar.LogLik;
        GradStar = ResStar.Grad;

        if Parameters.DiffusionType == 'Add'
            if Parameters.EstimatedVariables.SigmaRW
                qStarTemp = normpdf(SigmaRWStar, MuTemp, SigmaTemp);
                MuStar = SigmaRWStar + epsil^2*GStar^-1*GradStar/2;
                SigmaStar = epsil*sqrt(GStar^-1);
                qTempStar = normpdf(SigmaRWTemp, MuStar, SigmaStar);
            end
        end  
        LogRand = log(rand(1,1));
        LogAccRate = LogLikStar + log(qTempStar) - LogLikTemp - log(qStarTemp);

        if LogRand < LogAccRate
            if Parameters.DiffusionType == 'Add'
                if Parameters.EstimatedVariables.SigmaRW
                    ResTemp = ResStar;
                    LogLikTemp = LogLikStar;
                    GradTemp = GradStar;
                    SigmaRWTemp = SigmaRWStar;
                    GTemp = GStar;
                    MuTemp = MuStar;
                    SigmaTemp =SigmaStar;
                end
            end
            Accepted(IndIt) = 1;
        else
            Accepted(IndIt) = 0;
        end
        if Parameters.DiffusionType == 'Add'
            if Parameters.EstimatedVariables.SigmaRW
                Thetas(IndIt) = SigmaRWTemp;
            end
        end
        LogLiks(IndIt) = LogLikTemp;
        disp(['Acceptance Rate: ' num2str(sum(Accepted)/length(Accepted),2) ' (it' num2str(IndIt) ')'])
    end
    AccRate = sum(Accepted)/length(Accepted);
    if AccRate > 0.65
        epsil = epsil + 0.1;
    elseif AccRate < 0.55
        epsil = epsil - 0.1;
        if epsil<=0
            die
        end
    else
        Condition = 1;
    end
    disp(['For epsil = ' num2str(epsil) ', AccRate = ' num2str(AccRate)])
end


for IndIt = 1:1000
    if Parameters.DiffusionType == 'Add'
        if Parameters.EstimatedVariables.SigmaRW 
            SigmaRWStar = MuTemp + SigmaTemp*randn(1,1);
            ParametersStar = ParametersTemp;
            ParametersStar.SigmaRW = SigmaRWStar;
            GStar = (2*N)/SigmaRWStar^2;
        end
    end
    ResStar = EstimationSMCsmoothAndGrad(Data,ParametersStar);
    LogLikStar = ResStar.LogLik;
    GradStar = ResStar.Grad;
    
    if Parameters.DiffusionType == 'Add'
        if Parameters.EstimatedVariables.SigmaRW
            qStarTemp = normpdf(SigmaRWStar, MuTemp, SigmaTemp);
            MuStar = SigmaRWStar + epsil^2*GStar^-1*GradStar/2;
            SigmaStar = epsil*sqrt(GStar^-1);
            qTempStar = normpdf(SigmaRWTemp, MuStar, SigmaStar);
        end
    end  
    LogRand = log(rand(1,1));
    LogAccRate = LogLikStar + log(qTempStar) - LogLikTemp - log(qStarTemp);
    
    if LogRand < LogAccRate
        if Parameters.DiffusionType == 'Add'
            if Parameters.EstimatedVariables.SigmaRW
                ResTemp = ResStar;
                LogLikTemp = LogLikStar;
                GradTemp = GradStar;
                SigmaRWTemp = SigmaRWStar;
                GTemp = GStar;
                MuTemp = MuStar;
                SigmaTemp =SigmaStar;
            end
        end
        Accepted(IndIt) = 1;
    else
        Accepted(IndIt) = 0;
    end
    if Parameters.DiffusionType == 'Add'
        if Parameters.EstimatedVariables.SigmaRW
            Thetas(IndIt) = SigmaRWTemp;
        end
    end
    LogLiks(IndIt) = LogLikTemp;
    disp(['Acceptance Rate: ' num2str(sum(Accepted)/length(Accepted),2) ' (it' num2str(IndIt) ')'])
end


res.Data = Data;
res.Parameters = Parameters;
res.Thetas = Thetas;
res.LogLiks = LogLiks;
res.AccRate = sum(Accepted)/length(Accepted);
res.Espil = epsil;