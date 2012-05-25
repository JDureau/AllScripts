function [Parameters, TempPar] = CalibrateMethod( Data, Model, Parameters, TempPar,Mode)

if nargin == 4
    Mode = 'PMCMC';
end

Persistance = 0;
Epsils = [];
AccRates = [];

if strcmp(Parameters.MCMCType,'Lang')
    aim = 0.7;
    if (Parameters.ComputeRWsamplCov)
        if strcmp(Mode,'PMCMC') 
            Res = RunEstimationMethod(Data,Parameters,TempPar,Parameters.NbItRWCal);
        elseif strcmp(Mode,'MCMC') 
            Res = RunMCMC(Data, Model, Parameters, Parameters.NbItRWCal);
        end
        Parameters.RWsamplCov = cov(Res.LogThetas');
        Parameters.ComputeRWsamplCov = 0;
    end
elseif strcmp(Parameters.MCMCType,'Rand')
    aim = 0.23;
    if (Parameters.ComputeRWsamplCov)
         if strcmp(Mode,'PMCMC') 
            Res = RunEstimationMethod(Data,Parameters,TempPar,Parameters.NbItRWCal);
        elseif strcmp(Mode,'MCMC') 
            Res = RunMCMC(Data, Model, Parameters, Parameters.NbItRWCal);
         end
         Parameters.G = (cov(Res.LogThetas'))^-1;
        Parameters.ComputeRWsamplCov = 0;
    end
end
if not(isempty(Parameters.aim))
    aim = Parameters.aim;
end

        
NbItsRun = 100;
CoolingParameter = 0.98;
NbItsMax = 100;
Condition = 0;
it = 0;
Epsils = [];
AccRates = [];
EpsilsPot = [];
AccRatesPot = [];
TempAccRate = 1000;
EpsInit = 1;
while not(Condition)
    it = it+1;
    ParametersPot = Parameters;
    ParametersPot.Epsil = Parameters.Epsil + CoolingParameter^it*EpsInit*randn(1,1);
    if ParametersPot.Epsil > 0 
        if strcmp(Mode,'PMCMC') 
            Res = RunEstimationMethod(Data,Parameters,TempPar,NbItsRun);
        elseif strcmp(Mode,'MCMC') 
            Res = RunMCMC(Data, Model, Parameters, NbItsRun);
        end
        PotAccRate = Res.AccRate;
    else
        PotAccRate = 1000;
    end
    if abs(PotAccRate - aim)<abs(TempAccRate-aim)
        if strcmp(Mode,'PMCMC') 
            Res = RunEstimationMethod(Data,Parameters,TempPar,2*NbItsRun);
        elseif strcmp(Mode,'MCMC') 
            Res = RunMCMC(Data, Model, Parameters, 2*NbItsRun);
        end
        PotAccRate = Res.AccRate;
        if abs(PotAccRate - aim)<abs(TempAccRate-aim)
            TempAccRate = PotAccRate;
            Parameters = ParametersPot;
        end
        if abs(TempAccRate - aim) < 0.07
            Condition = 1;
        end
    end
    AccRates(end+1) = TempAccRate;
    Epsils(end+1) = Parameters.Epsil;
    AccRatesPot(end+1) = PotAccRate;
    EpsilsPot(end+1) = ParametersPot.Epsil;
    disp(['For epsil = ' num2str(Parameters.Epsil) ', AccRate = ' num2str(TempAccRate)])
    if  it >  NbItsMax
        Condition = 1;
    end
end


% Parameters.Calibrate{end+1} = [Epsils;AccRates];