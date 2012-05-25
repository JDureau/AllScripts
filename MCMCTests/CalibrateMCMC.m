function [Parameters, TempPar] = CalibrateMCMC( TempPar , Parameters)

Persistance = 0;
Epsils = [];
AccRates = [];

aim = Parameters.aim;
        
NbItsRun = 200;
CoolingParameter = 0.98;
NbItsMax = 500;
Condition = 0;
it = 0;
Epsils = [];
AccRates = [];
TempAccRate = 1000;
EpsInit = 1;
while not(Condition)
    it = it+1;
    ParametersPot = Parameters;
    ParametersPot.Epsil = Parameters.Epsil + CoolingParameter^it*EpsInit*randn(1,1);
    if ParametersPot.Epsil > 0 
        Res = RunMCMC(TempPar,ParametersPot,NbItsRun);
        PotAccRate = Res.AccRate;
    else
        PotAccRate = 1000;
    end
    if abs(PotAccRate - aim)<abs(TempAccRate-aim)
        Res = RunMCMC(TempPar,ParametersPot,2*NbItsRun);
        PotAccRate = Res.AccRate;
        if abs(PotAccRate - aim)<abs(TempAccRate-aim)
            TempAccRate = PotAccRate;
            Parameters = ParametersPot;
            disp(TempAccRate)
        end
        if abs(TempAccRate - aim) < 0.01
            Condition = 1;
        end
    end
    AccRates(end+1) = TempAccRate;
    Epsils(end+1) = Parameters.Epsil;
%     disp(['For epsil = ' num2str(Parameters.Epsil) ', AccRate = ' num2str(TempAccRate)])
    if  it >  NbItsMax
        disp('Didn''t converge')
        Condition = 1;
    end
end

Parameters.AccRate = TempAccRate;
disp(['Epsil: ' num2str(Parameters.Epsil)])
% Parameters.Calibrate{end+1} = [Epsils;AccRates];