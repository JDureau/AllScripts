function [] = CreateResGens()

temp = load('ParametersMysore.mat');
ParametersMysore = temp.Parameters;

temp = load('ParametersBelgaum.mat');
ParametersBelgaum = temp.Parameters;



ParametersBelgaum.NbVariables = 10;
ParametersBelgaum.SigmaObs = 0.1;
ParametersBelgaum.Problem = 'ImperialHIV';
ParametersBelgaum.DiffusionType = 'Affine';
ParametersBelgaum.ObservationLength = 25*12;
ParametersBelgaum.ComputationTStep = 0.5;

ParametersMysore.NbVariables = 10;
ParametersMysore.SigmaObs = 0.1;
ParametersMysore.Problem = 'ImperialHIV';
ParametersMysore.DiffusionType = 'Affine';
ParametersMysore.ObservationLength = 25*12;
ParametersMysore.ComputationTStep = 0.5;
Parameters = ParametersMysore;


% t0 = jan 85.
Data.Observations = zeros(10,5);
Data.Observations(7,2) = 26.11;
Data.Observations(7,3) = 24.24;
Data.Observations(8,4) = 5.4;
Data.Observations(7,5) = 11.10;
Data.Instants = round([0 236 264 286 292]/(Parameters.ComputationTStep));
Data.ObservedVariables = [ 0 7 7 8 7];
Data.NbComputingSteps = [0 diff(Data.Instants)];



HIVModel = struct();
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIVInitialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV_SMC_projection;


HIV2Model = struct();
HIV2Model.EKF_projection = @HIV2_EKF_projection;
HIV2Model.InitializeParameters = @HIV2Initialize;
HIV2Model.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
HIV2Model.SMC_projection = @HIV2_SMC_projection;

Parameters.TypeWork = 'Normal';

Parameters.ObsNoise = 0;
Ress = {};
Model = HIVModel;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: apply for a given set of sigmoids

xis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
Parameters.NbTSteps = length(xis);
Names = Parameters.Names.All;
for i = 1:length(Names)
    Parameters.(Names{i}).Sample = 1;
end




NbTests = 400;

ampls = [];
inflpts = [];
Ress = {};
cpt = 0;
ParsRecord = [];
ResGens1 = {};
ResGens2 = {};
ResGens3 = {};
ResGens4 = {};
while cpt < 800
    disp(cpt)
%    if 1
    try
        NbTests = 1000;       
        ampl = ((0.85*rand(1,1))^(0.8));
        baseline = max(0.1,rand(1,1)*(1-ampl));
        asympt = min(0.9,baseline+ampl);
        ampl = asympt-baseline;
        inflpt = rand(1,1)^(0.98)*450+100;
        steepness = rand(1,1)*30+20;
%
 %       Data.Instants = round([0 120 236 264 286 292]/(Parameters.ComputationTStep));
  %      Data.ObservedVariables = [ 0 7 7 7 8 7];
   %     Data.NbComputingSteps = [0 diff(Data.Instants)];
        
        xis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
        Parameters = SampleParameters(Parameters);
        Fts = baseline+ampl*Sigmoid((xis-inflpt*Parameters.ComputationTStep)/steepness);
        Parameters.MultNoise = 0.05;
        ResGen = GenerateDataForGivenSig(Data,Parameters,HIV2Model,Fts);
        if(Fts(end)-Fts(1)<0.05)
            if rand(1,1)<0.1
                die
            end
        end
%         if(inflpt<200)
%             if rand(1,1)<0.75
%                 die
%             end
%         end

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
	
  
	RestrictedParameters = DefineEstimatedParametersIndexes(RestrictedParameters);
	RestrictedParameters = DefineTransfFunctions(RestrictedParameters);
	RestrictedParameters = DefinePriors(RestrictedParameters);
	RestrictedParameters = UpdateParsNoTransfToTransf(RestrictedParameters);
%	ResGen = GenerateDataForGivenSig(Data,RestrictedParameters,HIV2Model,Fts);	


	ResGen.RestrictedParameters = RestrictedParameters;

        BiasedParameters = DefineEstimatedParametersIndexes(BiasedParameters);
        BiasedParameters = DefineTransfFunctions(BiasedParameters);
        BiasedParameters = DefinePriors(BiasedParameters);
        BiasedParameters = UpdateParsNoTransfToTransf(BiasedParameters);
	ResGen.BiasedParameters = BiasedParameters;
        
        cpt = cpt+1;
        if cpt<201
            ResGens1{end+1} = ResGen;
        elseif cpt<401
            ResGens2{end+1} = ResGen;
	elseif cpt<601
            ResGens3{end+1} = ResGen;
	else
            ResGens4{end+1} = ResGen;
        end
    end
end

save(['VIH_PlayingWithSigm_ResGens1.mat'],'ResGens1')
save(['VIH_PlayingWithSigm_ResGens2.mat'],'ResGens2')
save(['VIH_PlayingWithSigm_ResGens3.mat'],'ResGens3')
save(['VIH_PlayingWithSigm_ResGens4.mat'],'ResGens4')
%save(['VIH_PlayingWithSigm_ResGens1.mat'],'ResGens1')
%save(['VIH_PlayingWithSigm_ResGens2.mat'],'ResGens2')
%save(['VIH_PlayingWithSigm_ResGens3_Restricted.mat'],'ResGens3')
%save(['VIH_PlayingWithSigm_ResGens4_Restricted.mat'],'ResGens4')



