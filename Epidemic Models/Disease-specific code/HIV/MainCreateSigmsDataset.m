function [] = MainCreateSigmsDataset(ind)

cd('/users/ecologie/dureau/src/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])

SavePath = '/users/ecologie/dureau/src/AllData/Avahan';
load([SavePath '/ParametersMysore.mat']);


Parameters.Problem = 'ImperialHIV';
Parameters.SigmaObs = 0.1;
Parameters.ComputationTStep = 0.5;
Parameters.ObservationLength = 25*12;
Parameters.InitialFt.Sample = 1;
Parameters.InitialIPropF.Sample = 1;
Parameters.InitialIPropM.Sample = 1;
Parameters.TotMFactor.Sample = 1;
Parameters.Alpham1.Sample = 1;
Parameters.MuFm1.Sample = 1;
Parameters.MuMm1.Sample = 1;
Parameters.BetaMFPerAct.Sample = 1;
Parameters.BetaFMPerAct.Sample = 1;
Parameters.NumberActsPerClient.Sample =1;
Parameters.eHIV.Sample = 1;
Parameters.CF1.Sample = 1;
Parameters.CF2.Sample = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Data.Observations = zeros(10,5);
Data.Observations(7,2) = 26.11;
Data.Observations(7,3) = 24.24;
Data.Observations(8,4) = 5.4;
Data.Observations(7,5) = 11.10;
Data.Instants = round([0 236 264 286 292]/(Parameters.ComputationTStep));
Data.ObservedVariables = [ 0 7 7 8 7];
Data.NbComputingSteps = [0 diff(Data.Instants)];


HIVModel = struct();
HIVModel.EKF_projection = @HIV3_EKF_projection;
HIVModel.InitializeParameters = @HIV2Initialize;
HIVModel.SMC_projection = @HIV3_SMC_projection;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';


Mode = 2;
% 1 = Normal
% 2 = Biased
% 3 = Restricted
% 4 = AddingObs


NbTests = 400;

ampls = [];
inflpts = [];
Ress = {};
cpt = 0;
ParsRecord = [];
ResGens = {};
while cpt < 100
    disp(cpt)
    try
        
         NbTests = 1000;       
        ampl = ((0.85*rand(1,1))^(0.8));
        baseline = max(0.1,rand(1,1)*(1-ampl));
        asympt = min(0.9,baseline+ampl);
        ampl = asympt-baseline;
        inflpt = rand(1,1)^(0.98)*450+100;
        steepness = rand(1,1)*30+20;
        
        if Mode == 4
            Data.Instants = round([0 120 236 264 286 292]/(Parameters.ComputationTStep));
            Data.ObservedVariables = [ 0 7 7 7 8 7];
            Data.NbComputingSteps = [0 diff(Data.Instants)];
        end
        
        
        xis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
        Parameters = SampleParameters(Parameters);
        Fts = baseline+ampl*Sigmoid((xis-inflpt*Parameters.ComputationTStep)/steepness);
        Parameters.MultNoise = 0.05;
        ResGen = GenerateDataForGivenSig(Data,Parameters,HIVModel,Fts);
        if(Fts(end)-Fts(1)<0.05)
            if rand(1,1)<0.1
                die
            end
        end
        
        if Mode == 3
            Names = Parameters.Names.Estimated;
            RestrictedParameters = Parameters;
            for i = 1:length(Names)
                if strcmp(Parameters.(Names{i}).TransfType,'Logit')*(~(strcmp(Names{i},'SigmaRW')))*(~(strcmp(Names{i},'InitialFt')))*(~(strcmp(Names{i},'InitialIPropM')))*(~(strcmp(Names{i},'InitialIPropF')))*(~(strcmp(Names{i},'BetaFMPerAct')))*(~(strcmp(Names{i},'BetaMFPerAct')))
                    st = Parameters.(Names{i}).MaxLim-Parameters.(Names{i}).MinLim;
                    RestrictedParameters.(Names{i}).MinLim = RestrictedParameters.(Names{i}).MinLim+st/4;
                    RestrictedParameters.(Names{i}).MaxLim = RestrictedParameters.(Names{i}).MaxLim-st/4;
                    RestrictedParameters.(Names{i}).Value = (RestrictedParameters.(Names{i}).MaxLim+RestrictedParameters.(Names{i}).MinLim)/2;
                end
            end
            RestrictedParameters = DefineEstimatedParametersIndexes(RestrictedParameters);
            RestrictedParameters = DefineTransfFunctions(RestrictedParameters);
            RestrictedParameters = DefinePriors(RestrictedParameters);
            RestrictedParameters = UpdateParsNoTransfToTransf(RestrictedParameters);
            ResGen.RestrictedParameters = RestrictedParameters;
        	ResGen = GenerateDataForGivenSig(Data,RestrictedParameters,HIVModel,Fts);	

        elseif Mode == 2
            Names = Parameters.Names.Estimated;
            BiasedParameters = Parameters;
            for i = 1:length(Names)
                if strcmp(Parameters.(Names{i}).TransfType,'Logit')*not(strcmp(Names{i},'SigmaRW'))*not(strcmp(Names{i},'InitialFt'))*not(strcmp(Names{i},'eHIV'))
                    m = (Parameters.(Names{i}).MinLim+Parameters.(Names{i}).MaxLim)/2;
                    st = (Parameters.(Names{i}).MaxLim-Parameters.(Names{i}).MinLim);
                    if m>4/5*st
                        if rand(1,1)>0.5
                            swtch  = 1;
                        else
                            swtch  = -1;
                        end
                    else                  
                        swtch = 1;
                    end
                    BiasedParameters.(Names{i}).MinLim = max(eps,BiasedParameters.(Names{i}).MinLim + swtch*st/4);
                    BiasedParameters.(Names{i}).MaxLim = max(2*eps,BiasedParameters.(Names{i}).MaxLim + swtch*st/4);
                    newst = (BiasedParameters.(Names{i}).MaxLim-BiasedParameters.(Names{i}).MinLim);
                    if newst<0.9*st
                        disp(newst)
                        disp(st)
                        die
                    end
                    BiasedParameters.(Names{i}).Value = min(BiasedParameters.(Names{i}).MaxLim-0.1*newst,max(BiasedParameters.(Names{i}).MinLim+0.1*newst,BiasedParameters.(Names{i}).Value));
                end
            end
            BiasedParameters = DefineEstimatedParametersIndexes(BiasedParameters);
            BiasedParameters = DefineTransfFunctions(BiasedParameters);
            BiasedParameters = DefinePriors(BiasedParameters);
            BiasedParameters = UpdateParsNoTransfToTransf(BiasedParameters);
            ResGen.BiasedParameters = BiasedParameters;
        end
	
  
	


	

       
        
        
      
%         if(inflpt<200)
%             if rand(1,1)<0.75
%                 die
%             end
%         end
%         Names = Parameters.Names.Estimated;
%         for i = 1:length(Names)
%             ParsRecord(cpt+1,i) = Parameters.(Names{i}).Value;
%         end
%         TellParsValues(ResGen.Parameters)
%         ampls(end+1) = Fts(end)-Fts(1);
%         inflpts(end+1) = inflpt;
            cpt = cpt+1;
            ResGens{end+1} = ResGen;
%         pause(0.01)
    end
end
%subplot(2,1,1)
%hist(inflpts)
%subplot(2,1,2)
%plot(inflpts,ampls,'.')

Modes = {'','_Biased','_Restricted','_AddingObs'};


save([SavePath '/VIH_PlayingWithSigm_ResGens' Modes{Mode} '_' num2str(ind+1) '.mat'],'ResGens')

