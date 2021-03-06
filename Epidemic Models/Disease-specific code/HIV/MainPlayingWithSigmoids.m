% Main Sigm process HIV



%% load paths
%
cd('/Users/dureaujoseph/Dropbox/Taf/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])



addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Disease-specific code')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Disease-specific code\HIV')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Toolboxes\Resampling\pf_resampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Generic PMCMC tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Epidemic Models\Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts From Macbook\Toolboxes')



%%

% SavePath = '/Users/dureaujoseph/Documents/Taf/These/Results/Avahan';
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan';
temp = load([SavePath '/ParametersMysore.mat']);
ParametersMysore = temp.Parameters;


% SavePath = 'S:\Results';
temp = load([SavePath '/ParametersBelgaum.mat']);
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
ParametersMysore.DiffusionType = 'Logistic';
ParametersMysore.ObservationLength = 25*12;
ParametersMysore.ComputationTStep = 0.5;
Parameters = ParametersMysore;

% t0 = jan 85.
% Data.Observations = zeros(10,5);
% Data.Observations(7,2) = 26.11;
% Data.Observations(7,3) = 24.24;
% Data.Observations(8,4) = 5.4;
% Data.Observations(7,5) = 11.10;
% Data.Instants = round([0 50 100 150 200]/(Parameters.ComputationTStep));
% Data.ObservedVariables = [ 0 7 7 8 7];
% Data.NbComputingSteps = [0 diff(Data.Instants)];


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
HIVModel.EKF_projection = @HIV3_EKF_projection;
HIVModel.InitializeParameters = @HIVInitialize;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
HIVModel.SMC_projection = @HIV3_SMC_projection;


HIV2Model = struct();
HIV2Model.EKF_projection = @HIV2_EKF_projection;
HIV2Model.InitializeParameters = @HIV2Initialize;
HIV2Model.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';
HIV2Model.SMC_projection = @HIV2_SMC_projection;


%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Testing here the logistic regression
Parameters.CUinit.Value = 0.2;
Parameters.CUinit.Min = -10^14;
Parameters.CUinit.Max = 10^14;
Parameters.CUinit.MaxLim = 0.99;
Parameters.CUinit.MinLim = 0.01;
Parameters.CUinit.Estimated = 1;
Parameters.CUinit.TransfType = 'Logit';
Parameters.CUinit.Init = 1;

Parameters.CUdelta.Value = 0.2;
Parameters.CUdelta.Min = -10^14;
Parameters.CUdelta.Max = 10^14;
Parameters.CUdelta.MaxLim = 0.99;
Parameters.CUdelta.MinLim = 0.01;
Parameters.CUdelta.Estimated = 1;
Parameters.CUdelta.TransfType = 'Logit';
Parameters.CUdelta.Init = 1;

Parameters.k.Value = 0.2;
Parameters.k.Min = -10^14;
Parameters.k.Max = 10^14;
Parameters.k.Estimated = 1;
Parameters.k.TransfType = 'Log';
Parameters.k.Init = 1;


HIVapplyInference(ResGen.Data,ResGen.Parameters,Model);


Parameters.InitialFt.Value = 0.5;
Parameters.InitialFt.Min = -10^14;
Parameters.InitialFt.Max = 10^14;
Parameters.InitialFt.MaxLim = 0.93;
Parameters.InitialFt.MinLim = 0.07;
Parameters.InitialFt.Estimated = 1;
Parameters.InitialFt.TransfType = 'Logit';
Parameters.InitialFt.Init = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SavePath = 'S:\Results';
% Res = load([SavePath '\HIV_Mysore_PlayingPriors.mat']);
% Res = Res.Ress;
% Res = Res{17,3};


% Parameters = Res.Parameters;
% Names = Parameters.Names.All;
% for i = 1:length(Names)
%     Parameters.(Names{i}).Estimated = 1;
% end
% Parameters = DefineEstimatedParametersIndexes(Parameters);
% Parameters = DefineTransfFunctions(Parameters);
% Parameters = DefinePriors(Parameters);
% Parameters = UpdateParsNoTransfToTransf(Parameters);
Parameters.TypeWork = 'Normal';

Parameters.ObsNoise = 0;
Ress = {};
Model = HIVModel;
% 
% trajs = [];
% ampls = [];
% for indTest = 1:NbTests
%     disp(indTest)
%     ResGen = GenerateDataForGivenSig(Data,Parameters,HIVModel);
%     trajs(indTest,:) = ResGen.Data.BuiltTraj(:,9);    
%     ampls(indTest) = ResGen.Data.BuiltTraj(end,9)-ResGen.Data.BuiltTraj(1,9);
% end
% subplot(2,1,1)
% plot(trajs')
% subplot(2,1,2)
% hist(ampls)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step 1: apply for a given set of sigmoids

xis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
Parameters.NbTSteps = length(xis);
Names = Parameters.Names.All;
for i = 1:length(Names)
    Parameters.(Names{i}).Sample = 1;
end

% constant trajectories
NbTests = 20;
ampl = 0;
baselines = rand(1,NbTests);
asympts = baselines+ampl;
inflpt = 0;
steepness = 0;

Ress = {};
for i = 1 : NbTests
    i
    baseline = baselines(i);
    asympt = asympts(i);
    Parameters = SampleParameters(Parameters);
    Fts = baseline+ampl*Sigmoid((xis-inflpt)/steepness);
    try
        ResGen = GenerateDataForGivenSig(Data,Parameters,Model,Fts);
        Res = HIVapplyInference(ResGen.Data,ResGen.Parameters,Model);
        Ress{end+1} = Res;
        'OK'
    end
end


SavePath = 'S:\Results\';
save([SavePath 'VIH_PlayingWithSigm_Examples_Constants.mat'],'Ress')



% Early shifts
NbTests = 1000;
CptMax = 3;
ampls = rand(1,NbTests);
baselines = rand(1,NbTests).*(1-ampls);
asympts = baselines+ampls;
inflpt = 60;
steepnesses = rand(1,NbTests)*50;

Ress = {};
cpt = 0;
it = 0;
while and(cpt< 1,it<NbTests)
    it = it+1;
    cpt
    baseline = baselines(it);
    asympt = asympts(it);
    Parameters = SampleParameters(Parameters);
    Fts = baseline+ampls(it)*Sigmoid((xis-inflpt)/steepnesses(it));
    try
        Parameters.SwitchToIBM = 1;
        ResGen = GenerateDataForGivenSig(Data,Parameters,Model,Fts);
        Res = HIVapplyInference(ResGen.Data,ResGen.Parameters,Model);
        Ress{end+1} = Res;
        'OK'
        cpt = cpt+1;
    end
end


SavePath = 'S:\Results\';
save([SavePath 'VIH_PlayingWithSigm_Examples_EarlyShifts.mat'],'Ress')





% Middle shifts
NbTests = 1000;
ampls = rand(1,NbTests);
baselines = rand(1,NbTests).*(1-ampls);
asympts = baselines+ampls;
inflpt = 160;
steepnesses = rand(1,NbTests)*50;

Ress = {};
cpt = 0;
it = 0;
while and(cpt< 10,it<NbTests)
    it = it+1;
    cpt
    baseline = baselines(it);
    asympt = asympts(it);
    Parameters = SampleParameters(Parameters);
    Fts = baseline+ampls(it)*Sigmoid((xis-inflpt)/steepnesses(it));
    try
        ResGen = GenerateDataForGivenSig(Data,Parameters,Model,Fts);
        Res = HIVapplyInference(ResGen.Data,ResGen.Parameters,Model);
        Ress{end+1} = Res;
        'OK'
        cpt = cpt+1;
    end
end


SavePath = 'S:\Results\';
save([SavePath 'VIH_PlayingWithSigm_Examples_MiddleShifts.mat'],'Ress')




% Late shifts
NbTests = 1000;
ampls = rand(1,NbTests);
baselines = rand(1,NbTests).*(1-ampls);
asympts = baselines+ampls;
inflpt = 250;
steepnesses = rand(1,NbTests)*50;

Ress = {};
cpt = 0;
it = 0;
while and(cpt< 10,it<NbTests)
    it = it+1;
    cpt
    baseline = baselines(it);
    asympt = asympts(it);
    Parameters = SampleParameters(Parameters);
    Fts = baseline+ampls(it)*Sigmoid((xis-inflpt)/steepnesses(it));
    try
        ResGen = GenerateDataForGivenSig(Data,Parameters,Model,Fts);
        Res = HIVapplyInference(ResGen.Data,Parameters,Model);
        Ress{end+1} = Res;
        'OK'
        cpt = cpt+1;
    end
end


SavePath = 'S:\Results\';
save([SavePath 'VIH_PlayingWithSigm_Examples_LateShifts.mat'],'Ress')

SavePath = 'S:\Results\';
load([SavePath 'VIH_PlayingWithSigm_Examples_Constants.mat'])
load([SavePath 'VIH_PlayingWithSigm_Examples_EarlyShifts.mat'])
load([SavePath 'VIH_PlayingWithSigm_Examples_MiddleShifts.mat'])
load([SavePath 'VIH_PlayingWithSigm_Examples_LateShifts.mat'])

clf
for i = 1 : length(Ress)
    Res = Ress{i};
    Pars = Res.Parameters;
    Pars.TypeWork = 'Boston Examples';
    PlotResHIV(Res,Pars)
    i
    pause()
end



%% All kinds of tests


NbTests = 400;

ampls = [];
inflpts = [];
Ress = {};
cpt = 0;
ParsRecord = [];
ResGens1 = {};
ResGens2 = {};
while cpt < 50
    disp(cpt)
    try
        NbTests = 1000;       
        ampl = ((0.85*rand(1,1))^(0.8));
        baseline = max(0.1,rand(1,1)*(1-ampl));
        asympt = min(0.9,baseline+ampl);
        ampl = asympt-baseline;
        inflpt = rand(1,1)^(0.98)*450+100;
        steepness = rand(1,1)*30+20;

%          Data.Instants = round([0 120 236 264 286 292]/(Parameters.ComputationTStep));
%         Data.ObservedVariables = [ 0 7 7 7 8 7];
%         Data.NbComputingSteps = [0 diff(Data.Instants)];
        
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
        for i = 1:length(Names)
            ParsRecord(cpt+1,i) = Parameters.(Names{i}).Value;
        end
        TellParsValues(ResGen.Parameters)
        ampls(end+1) = Fts(end)-Fts(1);
        inflpts(end+1) = inflpt;
        cpt = cpt+1;
        if cpt<201
            ResGens1{end+1} = ResGen;
        else
            ResGens2{end+1} = ResGen;
        end
        pause(0.01)
    end
end

save([SavePath 'VIH_PlayingWithSigm_ResGens1.mat'],'ResGens1')
save([SavePath 'VIH_PlayingWithSigm_ResGens2.mat'],'ResGens2')

inds1 = find(inflpts<240);
inds2 = find(and(240<inflpts,inflpts<456));
inds3 = find(456<inflpts);

disp(num2str([length(inflpts) length(unique(inflpts))]))
disp(num2str([length(inds1) length(inds2) length(inds3)]))

subplot(3,1,1)
hist(ampls,100)
subplot(3,1,2)
hist(inflpts)
subplot(3,1,3)
plot(ampls,inflpts,'.')


k = ceil(sqrt(length(Names)));
for i = 1:length(Names)
    subplot(k,k,i)
    plot(ParsRecord(:,i),inflpts,'.')
%     [fi,xi] = ksdensity(ParsRecord(inds1,i));
%     plot(xi,fi,'r')
%     hold on
%     [fi,xi] = ksdensity(ParsRecord(inds2,i));
%     plot(xi,fi,'g')
%     [fi,xi] = ksdensity(ParsRecord(inds3,i));
%     plot(xi,fi,'b')
%     hold off
end

SavePath = 'S:\Results\';
load([SavePath 'VIH_PlayingWithSigm_ResGens.mat'])

inds = 26:50;
Ress = {};
for i = inds(1):inds(end)
%     Names = Parameters.Names.Estimated;
%     for i = 1:length(Names)
%        ParametersBelgaum.(Names{i}).Value = min(max(ResGen.Parameters.(Names{i}).Value,ParametersBelgaum.(Names{i}).MinLim*1.02),ParametersBelgaum.(Names{i}).MaxLim*0.98);
%        ParametersBelgaum.(Names{i}).Estimated = 1;
%     end
%     ParametersBelgaum = UpdateParsNoTransfToTransf(ParametersBelgaum);
%     ParametersBelgaum.MultNoise = Parameters.MultNoise;
             
    Res = HIVapplyInference(ResGens{i}.Data,ResGens{i}.Parameters);
%         Res = HIVapplyInference(ResGen.Data,ParametersBelgaum);
% 
    Ress{end+1} = Res;
    %         
    %         
    SavePath = 'S:\Results\';
    save([SavePath 'VIH_PlayingWithSigm_RightPriors2.mat'],'Ress')
end


ampls = [];
inflpts = [];
Ress = {};
cpt = 0;
ParsRecord = [];
while cpt < 40
    disp(cpt)
    try
        NbTests = 1000;       
        ampl = ((0.85*rand(1,1))^(0.8));
        baseline = max(0.1,rand(1,1)*(1-ampl));
        asympt = min(0.9,baseline+ampl);
        ampl = asympt-baseline;
        inflpt = rand(1,1)^(0.9)*400+150;
        steepness = rand(1,1)*30+20;

        xis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
        Parameters = SampleParameters(Parameters);
        Fts = baseline+ampl*Sigmoid((xis-inflpt*Parameters.ComputationTStep)/steepness);
        
%         Data.Instants = round([0 120 236 264 286 292]/(Parameters.ComputationTStep));
%         Data.ObservedVariables = [ 0 7 7 7 8 7];
%         Data.NbComputingSteps = [0 diff(Data.Instants)];

        
%         Parameters.SwitchToIBM = 1;
        Parameters.MultNoise = 0.05;
        ResGen = GenerateDataForGivenSig(Data,Parameters,HIV2Model,Fts);
        if(Fts(end)-Fts(1)<0.05)
            if rand(1,1)<0.1
                die
            end
        end
        if(inflpt<200)
            if rand(1,1)<0.75
                die
            end
        end
        TellParsValues(ResGen.Parameters)
        ampls(end+1) = Fts(end)-Fts(1);
        inflpts(end+1) = inflpt;
        
%         Names = Parameters.Names.Estimated;
%         for i = 1:length(Names)
%             ParametersBelgaum.(Names{i}).Value = min(max(ResGen.Parameters.(Names{i}).Value,ParametersBelgaum.(Names{i}).MinLim*1.02),ParametersBelgaum.(Names{i}).MaxLim*0.98);
%             ParametersBelgaum.(Names{i}).Estimated = 1;
%         end
%         ParametersBelgaum = UpdateParsNoTransfToTransf(ParametersBelgaum);
%         ParametersBelgaum.MultNoise = Parameters.MultNoise;
%         ResGen.Parameters.SwitchToIBM = 0;
% %         
%         Res = HIVapplyInference(ResGen.Data,ResGen.Parameters);
% %         Res = HIVapplyInference(ResGen.Data,ParametersBelgaum);
% % 
%         Ress{end+1} = Res;
% %         
% %         
%         SavePath = 'S:\Results\';
%         save([SavePath 'VIH_PlayingWithSigm_RightPriors.mat'],'Ress')
        cpt = cpt+1;
    end
end

subplot(3,1,1)
hist(ampls,100)
subplot(3,1,2)
hist(inflpts)
subplot(3,1,3)
plot(ampls,inflpts,'.')

Resss = {};
inds = {'','2','3','4','5','6'};
for i = 1:6%3
    SavePath = 'S:\Results\';
    load([SavePath 'VIH_PlayingWithSigm_RightPriors' inds{i} '.mat'])
    for j = 1:length(Ress)
        Resss{end+1} = Ress{j};
        length(Resss)
    end
end


MeanErrors = [];
Errors = [];
amplin = [];
amplout = [];
slope = [];
inflptin = [];
inflptout = [];
corrs = [];
for j = 1:length(Resss)
    j
    Res = Resss{j};
    if and(Res.AccRate > 0.18, Res.AccRate < 0.40)
        tmp = diff(Res.Data.Fts);
        [b,ind] = max(tmp);
        if b<0.1*Parameters.ComputationTStep
            tmp = squeeze(Res.Paths(:,3,:));
            tmp = exp(tmp)./(1+exp(tmp));
            if 1%and(ind>50,ind<530)
                Pars = FitSigmoid(mean(tmp),Res.Parameters);
                amplin(end+1) = Res.Data.Fts(585)-Res.Data.Fts(1)   ;
%                 plot(Pars.Sigm,'k')
%                 hold on
%                 plot(mean(tmp),'g')
%                 hold off     
%                 pause(0.01)
                amplout(end+1) = Pars.Ampl;           
                inflptin(end+1) = ind;
                [b,indmed] = min(abs(Pars.Sigm-(Pars.Sigm(585)+Pars.Sigm(1))/2));
                inflptout(end+1) = indmed;
%                 if length(inflptout) == 63
                    Res.Parameters.TypeWork='Boston Examples HIV2';
                    Res.Parameters.Sigm = Pars.Sigm; 
                    PlotResHIV(Res,Res.Parameters,Res.Data.Fts)
                    title(num2str(length(inflptout)))
%                   
%                     die
%                 end
                pause(0.01)
            end
        end
    end
end
            
            
IndsOut= [1 4 20 24 26 21 52 63 64 73];
inds2(IndsOut) = [];

Res = Resss{4};
Pars = Res.Parameters;
Pars.TypeWork = 'Boston Examples';
PlotResHIV(Res,Pars)

list = [60 73 4 44];
BetterRess = {};
for i = 1:length(list)
    BetterRess{i} = HIVapplyInference(Resss{list(i)}.Data,Resss{list(i)}.Parameters,Model); 
end

AccRates = [];
for j = 1:length(Resss)
    if Resss{j}.AccRate>0.15
        Res = Resss{j};
        Names = Res.Parameters.Names.Estimated;
        for i = 1:length(Names)
            subplot(ceil(sqrt(length(Names))),ceil(sqrt(length(Names))),i)
            plot(Res.TransfThetas(i,:))
            title(Names{i})
        end
        pause()
    end
end


SavePath = 'S:\Results\';
save([SavePath 'VIH_PlayingWithSigm_4exsWithMoreIts.mat'],'BetterRess')

Res = BetterRess{3};
Pars = Res.Parameters;
Pars.TypeWork = 'Boston Examples';
PlotResHIV(Res,Pars)


plot(amplin,MeanErrors,'.')
inds1 = find(amplin>0.4);
inds2 = find(amplin<=0.4);
plot(amplin(inds1),MeanErrors(inds1),'.r')
hold on
plot(amplin(inds2),MeanErrors(inds2),'.g')
hold off

% inflt inflt
clf
plot(0:0.01:600,0:0.01:600,'--k','LineWidth',4)
hold on
for j = 1:length(inflptin)
    if inflptin(j) < 240
        if amplin(j)<0.4
            scatter(inflptin(j),inflptout(j),25,'r','filled','MarkerEdgeColor','k')
        else
            scatter(inflptin(j),inflptout(j),95,'r','filled','MarkerEdgeColor','k')
        end
    elseif inflptin(j) < 432
        if amplin(j)<0.4
            scatter(inflptin(j),inflptout(j),25,'g','filled','MarkerEdgeColor','k')
        else
            scatter(inflptin(j),inflptout(j),95,'g','filled','MarkerEdgeColor','k')
        end
    elseif inflptin(j) < 600
        if amplin(j)<0.4
            scatter(inflptin(j),inflptout(j),25,'b','filled','MarkerEdgeColor','k')
        else
            scatter(inflptin(j),inflptout(j),95,'b','filled','MarkerEdgeColor','k')
        end
    end
end
% scatter(inflptin(inds1),inflptout(inds1),95,'g','filled','MarkerEdgeColor','k')
% scatter(inflptin(inds2),inflptout(inds2),25,'g','filled','MarkerEdgeColor','k')
for j = 1:length(inflptin)
    text(inflptin(j),inflptout(j),num2str(j))
end
hold off
xlabel('Simulated CU inflection point','Fontsize',16)
ylabel('Estimated CU inflection point','Fontsize',16)
xlim([0 600] )
ylim([0 600] )
title('Can we capture the inflection point?','FontWeight','bold','Fontsize',18)
set(gca,'XTick',[0:(600)/5:(600)])
set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
set(gca,'YTick',[0:(600)/5:(600)])
set(gca,'YTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])
corr(inflptin(inds2)',inflptout(inds2)')

% ampl ampl
Myscattercloud(amplin(inds2),amplout(inds2),20,0,1,0,1)
xlim([0 1])
ylim([0 1])
hold on
plot(0:0.01:1,0:0.01:1,'--','LineWidth',2)
hold off
xlabel('Input Amplitude','Fontsize',16)
ylabel('Output Amplitude','Fontsize',16)
title('How well is the imput amplitude captured?','FontWeight','bold','Fontsize',18)

clf
plot(0:0.01:1,0:0.01:1,'--k','LineWidth',4)
hold on
for j = 1:length(inflptin)
%     if and(amplin(j)<0.4,amplout(j)<0.4)
%         scatter(amplin(j),amplout(j),25,'g','filled','MarkerEdgeColor','k')
%     elseif and(amplin(j)>0.4,amplout(j)>0.4)
%         scatter(amplin(j),amplout(j),95,'g','filled','MarkerEdgeColor','k')
%     elseif and(amplin(j)>0.4,amplout(j)<0.4)
%         scatter(amplin(j),amplout(j),95,'r','filled','MarkerEdgeColor','k')
%     elseif and(amplout(j)>0.4,amplin(j)<0.4)
%         scatter(amplin(j),amplout(j),25,'r','filled','MarkerEdgeColor','k')
%     end
    if inflptin(j) < 240
        if amplin(j)<0.4
            scatter(amplin(j),amplout(j),25,'r','filled','MarkerEdgeColor','k')
        else
            scatter(amplin(j),amplout(j),95,'r','filled','MarkerEdgeColor','k')
        end
    elseif inflptin(j) < 432
        if amplin(j)<0.4
            scatter(amplin(j),amplout(j),25,'g','filled','MarkerEdgeColor','k')
        else
            scatter(amplin(j),amplout(j),95,'g','filled','MarkerEdgeColor','k')
        end
    elseif inflptin(j) < 600
        if amplin(j)<0.4
            scatter(amplin(j),amplout(j),25,'b','filled','MarkerEdgeColor','k')
        else
            scatter(amplin(j),amplout(j),95,'b','filled','MarkerEdgeColor','k')
        end
    end
end
% scatter(amplin(inds1),amplout(inds1),95,'g','filled','MarkerEdgeColor','k')
% scatter(amplin(inds2),amplout(inds2),25,'g','filled','MarkerEdgeColor','k')
% for j = 1:length(inflptin)
%     text(amplin(j)*(1+0.2*rand(1,1)),amplout(j),num2str(j))
% end
colormap('default')
plot(0.4*ones(2,1),[0 1],'--k','LineWidth',2)
plot([0 1],0.4*ones(2,1),'--k','LineWidth',2)
hold off
xlabel('Simulated CU amplitude','Fontsize',16)
ylabel('Estimated CU amplitude','Fontsize',16)
title('Can we capture the amplitude of the change in CU?','FontWeight','bold','Fontsize',18)
xlim([0 1] )
ylim([0 1] )
corr(amplin',amplout')

% baseline baseline
inds = find(and(baselineout(inds2)<1,0<baselineout(inds2)))
plot(0:0.01:1,0:0.01:1,'--k','LineWidth',2)
hold on
scatter(baselinein(tmpinds),baselineout(tmpinds),95,'g','filled','MarkerEdgeColor','k')
% scatter(baselinein(inds2(inds)),baselineout(inds2(inds)),95,'g','filled','MarkerEdgeColor','k')
hold off
xlabel('Input Baseline','Fontsize',16)
ylabel('Output Baseline','Fontsize',16)
% xlims([360:600] )
title('How well is the imput baseline captured?','FontWeight','bold','Fontsize',18)


[fi,xi] = max(baselinein(inds2(inds))./baselineout(inds2(inds)))
tmpinds = inds2(inds);
tmpinds(34) = [];

clf
[fi,xi] = ksdensity(MeanErrors(1,:));
plot(xi,fi)
hold on
[fi,xi] = ksdensity(MeanErrors(2,:));
plot(xi,fi,'k')
[fi,xi] = ksdensity(MeanErrors(3,:));
plot(xi,fi,'g')
hold off

[b,ind] = min(MeanErrors(3,:));
PlotResHIV(Resss{3}{ind},Resss{3}{ind}.Parameters,Resss{3}{ind}.Data.BuiltTraj(:,9))
[b,ind] = max(MeanErrors(3,:));
PlotResHIV(Resss{3}{ind},Resss{3}{ind}.Parameters,Resss{3}{ind}.Data.BuiltTraj(:,9))

hist([MeanErrors(3,:) MeanErrors(2,:)])
xlabel('Mean Squared error')
title('Distribution of Mean Squared Errors')

xs = [];
ys = [];
for j = 1:length(Resss{3})
    Res = Resss{i}{j};
    for k = 1:584
        xs(end+1) = k;
        ys(end+1) = Errors(3,j,k);
    end
end
scattercloud(xs*Parameters.ComputationTStep,ys)
xlabel('time')
ylabel('squared errors at time t')
title('Distribution of squared errors along time')

% role of ampl
plot(ampl(3,:),MeanErrors(3,:),'.')
hold on
plot(ampl(2,:),MeanErrors(2,:),'.')
plot(ampl(1,:),MeanErrors(1,:),'.')
hold off
xlabel('Amplitude of Ft trajectory')
ylabel('Mean Squarred Error')
title('Role of amplitude')

% role of t
plot(inflpt(3,:),MeanErrors(3,:),'.')
hold on
plot(inflpt(2,:),MeanErrors(2,:),'.')
plot(inflpt(1,:),MeanErrors(1,:),'.')
hold off
xlabel('Inflection point trajectory')
ylabel('Mean Squarred Error')
title('Role of inflection point')


PlotScatteredQuantity([([inflpt(3,:) inflpt(2,:) inflpt(1,:)]*Parameters.ComputationTStep)' [ampl(3,:)  ampl(2,:) ampl(1,:)]'],[corrs(3,:) corrs(2,:) corrs(1,:)],7)
% PlotScatteredQuantity([([inflpt(3,:) inflpt(2,:) inflpt(1,:)]*Parameters.ComputationTStep)' [ampl(3,:)  ampl(2,:) ampl(1,:)]'],-[MeanErrors(3,:) MeanErrors(2,:) MeanErrors(1,:)],6)
xlabel('Timing of inflection point in condom use')
ylabel('Amplitude of the change in condom use (in %)')
title('Combined role of amplitude and time of inflection on estimation quality (correlation)', 'FontWeight','bold')
colorbar
colorbar('YTickLabel', {'Worse','','','Average','','','Best',''},'FontWeight','bold')
set(gca,'XTick',[0:(50)/5:(50)])
set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])




plot([inflpt(3,:) inflpt(2,:) inflpt(1,:)],[ampl(3,:) ampl(2,:) ampl(1,:)],'.')

ind = 59;
traj = Resss{3}{ind}.Data.BuiltTraj(:,9);
RealDists = [];
EstDists  = [];
for i = 2:3
    for j = 1:100
        RealDists(end+1) = mean((traj-Resss{i}{j}.Data.BuiltTraj(:,9)).^2);
        errors = zeros(1,584);
        Res = Resss{i}{j};
        for k = 1:584
            errors(k) = mean((Res.Paths(:,3,k)-traj(k)).^2);
        end
        EstDists(end+1) = mean(errors);
    end
end

clf
plot(RealDists,EstDists,'.')
xlabel('Real Mean Squarred Distance to Avahan-like Trajectories')
ylabel('Estimated Mean Squarred Distance to Avahan-like Trajectories')
title('How well does this procedure discriminate trajectories like the ones estimated for Avahan?', 'FontWeight','bold')


RealInflPts = [];
RealEnd = [];
RealBaseLine = [];
EstInflPts = [];
EstEnd = [];
EstBaseLine = [];
InteresInds = [];
for i = 2:3
    for j = 1:100
        Res = Resss{i}{j};
        MeanTraj = mean(Res.Paths(:,3,:));
        RealInflPts(end+1) = inflpt(i,j);
        tmp = diff(MeanTraj);
        [b,ind] = max(tmp);
        EstInflPts(end+1) = ind;
        RealBaseLine(end+1) = Res.Data.BuiltTraj(1,9);
        EstBaseLine(end+1) = MeanTraj(1);
        RealEnd(end+1) = Res.Data.BuiltTraj(end,9);
        EstEnd(end+1) = MeanTraj(end);
        if and(ampl(i,j)>0.5,inflpt>420)
            InteresInds = length(EstEnd);
        end
    end
end


plot(RealDists,((-RealInflPts+EstInflPts)*Parameters.ComputationTStep),'.')
xlabel('Mean Squarred Distance to Avahan Estimated Trajectories')
ylabel('Error in inflection point estimation (in months) (Estimated-Real)')
title('How well is the inflection point estimated for Avahan-like trajectories?', 'FontWeight','bold')

plot(RealDists,abs((RealBaseLine-EstBaseLine)),'.')
xlabel('Mean Squarred Distance to Avahan Estimated Trajectories')
ylabel('Error in condom use baseline estimation ')
title('How well is the condom use baseline estimated for Avahan-like trajectories?', 'FontWeight','bold')

plot(RealDists,abs((RealEnd-EstEnd)),'.')
xlabel('Mean Squarred Distance to Avahan Estimated Trajectories')
ylabel('Error in condom use 2010 estimation ')
title('How well is the condom use in 2010 estimated for Avahan-like trajectories?', 'FontWeight','bold')





SavePath = 'S:\Results\';
load([SavePath 'VIH_PlayingWithSigm2.mat'])

MeanVals = [];
AmplCI = [];
Ress  = Resss{3};
for i = 1:3
    for j = 1:length(Resss{i})
        if ampl(i,j)<1
            Resss{i}{j}.Parameters.TypeWork = 'PlayObs';
            PlotResHIV(Resss{i}{j},Resss{i}{j}.Parameters,Resss{i}{j}.Data.BuiltTraj(:,9))
            title([i j])
            MeanVals(end+1) = mean(Resss{i}{j}.Data.BuiltTraj(:,9));
            AmplCI(end+1) = mean(quantile(Resss{i}{j}.Paths(:,3,:),0.75)-quantile(Resss{i}{j}.Paths(:,3,:),0.25));
            pause()
        end
    end        
end

i = 2
j = 26
PlotResHIV(Resss{i}{j},Resss{i}{j}.Parameters,Resss{i}{j}.Data.BuiltTraj(:,9))

clf
plot(abs(MeanVals-0.5),AmplCI,'.')
xlabel('Distance between condom use contant value and 0.5')
ylabel('Width of the mean posterior 50% confidence interval')
title('Estimation of constant condom use improve when this value is far from 0.5', 'FontWeight','bold')


hist(max(squeeze(Resss{1}{7}.Paths(:,3,:))')-min(squeeze(Resss{1}{7}.Paths(:,3,:))'))

PlotResHIV(Resss{1}{7},Resss{1}{7}.Parameters,Resss{1}{7}.Data.BuiltTraj(:,9))

hist(max(squeeze(Resss{2}{17}.Paths(:,3,:))')-min(squeeze(Resss{2}{17}.Paths(:,3,:))'))

PlotResHIV(Resss{2}{17},Resss{2}{17}.Parameters,Resss{2}{17}.Data.BuiltTraj(:,9))

clf
inds = ceil(rand(1,2)*100);
plot(squeeze(Resss{2}{17}.Paths(inds,3,:))','.')
title('fig 12: 2 random sample paths')

Parameters = Res.Parameters;
Trajs = {};

% cst low
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 1200;
Fts = 0.1*ones(1,Parameters.NbTSteps);
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 1;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;

%cst high: extinction
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 3000;
Fts = 0.5*ones(1,Parameters.NbTSteps);
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 2;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;

% sigmoid: realistic
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 100;
tmp = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ComputationTStep*Parameters.NbTSteps;
Fts = sigmoid((tmp-40)/2)*0.6 +0.01;
plot( Fts)
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 3;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;


% sigmoid: higher and slower
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 120;
tmp = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ComputationTStep*Parameters.NbTSteps;
Fts = sigmoid((tmp-40)/5)*0.6 +0.3;
plot( Fts)
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 4;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;

% random walk 
Parameters.ComputationTStep = 0.5;
Parameters.NbTSteps = 1200;
tmp = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ComputationTStep*Parameters.NbTSteps;
Fts = min(1,max(0,0.1+cumsum(0.01*randn(1,Parameters.NbTSteps)*sqrt(Parameters.ComputationTStep))));
plot( Fts)
HIV_CreateData(Fts,Parameters,HIVModel,Data)

ind = 5;
Trajs{ind} = struct();
Trajs{ind}.Fts = Fts;
Trajs{ind}.Parameters = Parameters;


ObsPatterns{1}.NbObs = 4;
ObsPatterns{1}.NbPatterns = 15;
ObsPatterns{2}.NbObs = 5;
ObsPatterns{2}.NbPatterns = 15;
for i = 1:2
    for j = 1:ObsPatterns{i}.NbPatterns 
        rdinds = sort(rand(1,ObsPatterns{i}.NbObs));
        rdclientorfsw = round(rand(1,ObsPatterns{i}.NbObs));
        for indtraj = 1:length(Trajs)
            Parameters = Trajs{indtraj}.Parameters;
            Fts = Trajs{indtraj}.Fts;
            ObsInds = ceil(rdinds*Parameters.NbTSteps);
            Data = HIV_CreateData(Fts,Parameters,HIVModel,Data);
            Data.Observations = zeros(10,ObsPatterns{i}.NbObs+2);
            Data.Instants = zeros(1,ObsPatterns{i}.NbObs+2);
            for indobs = 1:ObsPatterns{i}.NbObs
                Data.Observations(7+rdclientorfsw(indobs),indobs+1) = Data.BuiltTraj(ObsInds(indobs),7+rdclientorfsw(indobs));
                Data.Instants(indobs+1) = ObsInds(indobs);
                Data.ObservedVariables(indobs+1) = (7+rdclientorfsw(indobs));
            end
            Data.Observations(8,ObsPatterns{i}.NbObs+2)= 0;
            Data.Instants(ObsPatterns{i}.NbObs+2) = size(Data.BuiltTraj,1);
            Data.ObservedVariables(ObsPatterns{i}.NbObs+2) = 8; 
            Data.NbComputingSteps = [0 diff(Data.Instants)];
            Trajs{indtraj}.Datas{i,j} = Data;
            Model = HIVModel;
            tmp1 = zeros(1,10);
            tmp2 = zeros(1,10);
            tmp1(1,7) = 1;
            tmp2(1,8) = 1;
            temp = {};
            temp{1} = tmp1;
            temp{2} = tmp2;
            
            Model.ObservationJacobian = {};
            for indobs = 1:ObsPatterns{i}.NbObs
                Model.ObservationJacobian{1+indobs} = temp{1+rdclientorfsw(indobs)};
                Model.ObservationMeasurementNoise{1+indobs} = (Parameters.SigmaObs*Data.Observations(7+rdclientorfsw(indobs),1+indobs))^2;
            end
            Trajs{indtraj}.Models{i,j} = Model;
        end
    end
end
            


% go            
for indObsPattern = 1:2
    for indtraj = 1:length(Trajs)
        for j = 1:ObsPatterns{indObsPattern}.NbPatterns      
            Data = Trajs{indtraj}.Datas{indObsPattern,j};
            Model = Trajs{indtraj}.Models{indObsPattern,j};
            Parameters = Trajs{indtraj}.Parameters;
            Names = Parameters.Names.All;
            for i = 1:length(Names);
                Parameters.(Names{i}).Estimated = 0;
            end
            Parameters.InitialFt.Value = Trajs{indtraj}.Fts(1);
            Parameters.SigmaRW.Value = 0.1;
            Parameters.SigmaRW.Min = 0;
            Parameters.SigmaRW.Max = 10^14;
            Parameters.SigmaRW.Estimated = 1;
            Parameters.SigmaRW.TransfType = 'Log';
            Parameters = DefineEstimatedParametersIndexes(Parameters);
            Parameters = DefineTransfFunctions(Parameters);
            Parameters = UpdateParsNoTransfToTransf(Parameters);
            Parameters = DefinePriors(Parameters);
            % optimize RW parameters with SMC
            Parameters.NoPaths = 0;
            Names = Parameters.Names.Estimated;
            Initialization = [];
            for i = 1:length(Names)
                Initialization(i) = Parameters.(Names{i}).TransfValue ;
            end
            [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizewithPrior(x,Data,Model,Parameters),Initialization,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));
            Names = Parameters.Names.Estimated;
            for i = 1:length(Names)
                Parameters.(Names{i}).TransfValue = (x(i));
            end
            Parameters = UpdateParsTransfToNoTransf(Parameters);
            TellParsValues(Parameters)

            Parameters.TypeWork = 'PlayObs';
            Parameters.NbParticules = 5000;
            Parameters.DiffusionType = 'Add';
            Res = EstimationSMCfiltGen(Data, Model, Parameters);

            NbIts = 40;
            Paths = zeros(NbIts,length(Parameters.PathsToKeep),sum(Data.NbComputingSteps));
            LogLiks = [];
            for i = 1:NbIts
                disp(i)
                Temp = EstimationSMCsmoothGen(Data, Model, Parameters);
                ind = ceil(rand(1,1)*Parameters.NbParticules);
                Paths(i,:,:) = Temp.CompletePaths(ind,:,:);
                LogLiks(i) = Temp.LogLik;
            end
            Res.Paths = Paths;
            PlotResHIV(Res,Parameters)

            Trajs{indtraj}.Ress{indObsPattern,j} = Res;      
        end
    end
end


SavePath = 'S:\Results';
% SavePath = 'H:\PhD Work\Matlab Scripts\Epidemic Models\SIRS\Results';
% save([SavePath '\HIV_Mysore_ForProj_02I0_Calib.mat'],'ResRW1_RW_Calib_Rds_1_2_3')
save([SavePath '\HIV_PlayingObsPatterns.mat'],'Trajs')

load([SavePath '\HIV_PlayingObsPatterns.mat'])

clf
indtraj = 5;
indObsPattern = 2;
NbPats = 4;%ObsPatterns{indObsPattern}.NbPatterns ; 
tmp = ceil(NbPats/3);
Data = Trajs{indtraj}.Datas{indObsPattern,1};
TStep = Trajs{indtraj}.Parameters.ComputationTStep;
NbTSteps = Trajs{indtraj}.Parameters.NbTSteps;
xis = TStep:Parameters.ComputationTStep:NbTSteps*TStep;
subplot(tmp*6,1,1:tmp)
plot(xis,Data.BuiltTraj(:,7),'m')
title('FSWs')
subplot(tmp*6,1,tmp+1:2*tmp)
plot(xis,Data.BuiltTraj(:,8),'k')
title('Clients')
subplot(tmp*6,1,2*tmp+1:3*tmp)
plot(xis,Data.BuiltTraj(:,9))
ylim([0 1])
title('Condom Use')
scores = [];
for i = 1:ObsPatterns{indObsPattern}.NbPatterns
    Res = Trajs{indtraj}.Ress{indObsPattern,i};
    x = 0;
    for j = 1:size(Res.Paths,1)
        x = x + (squeeze(Res.Paths(j,3,:)) - Data.BuiltTraj(1:size(Res.Paths,3),9)).^2;
    end        
    scores(i) = mean(x);
end
[bof,inds] = sort(scores);
inds = inds([1 2 ObsPatterns{indObsPattern}.NbPatterns-1 ObsPatterns{indObsPattern}.NbPatterns])
for i = 1:4
    ind = inds(i);
    Data = Trajs{indtraj}.Datas{indObsPattern,ind};
    subplot(tmp*6,1,3*tmp+i)
    tmpinds = find(Data.ObservedVariables == 7);
    if not(isempty(tmpinds))
        plot(Data.Instants(tmpinds),0,'o','MarkerEdgeColor','m','MarkerFaceColor','m')
    end
    hold on
    tmpinds = find(Data.ObservedVariables == 8);
    if indObsPattern == 2
        tmpinds = tmpinds(1:end-1);
    else
        tmpinds = tmpinds(1:end-2);
    end
    if not(isempty(tmpinds))
        plot(Data.Instants(tmpinds),0,'o','MarkerEdgeColor','k','MarkerFaceColor','k')
    end
    hold off
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',0)
    set(gca,'yticklabel',scores(ind))
%     clf
%     plot(Data.BuiltTraj(1:size(Res.Paths,3),9),'g')
%     hold on
%     Res = Trajs{indtraj}.Ress{indObsPattern,ind};
%     plot(mean(squeeze(Res.Paths(:,3,:))))
%     hold off
%     pause()
end    


inds = [1 3 4 5];
for i = 1:4
    indtraj = inds(i);
    indObsPattern = 2;
    NbPats = 4;%ObsPatterns{indObsPattern}.NbPatterns ; 
    tmp = ceil(NbPats/3);
    Data = Trajs{indtraj}.Datas{indObsPattern,1};
    TStep = Trajs{indtraj}.Parameters.ComputationTStep;
    NbTSteps = Trajs{indtraj}.Parameters.NbTSteps;
    xis = TStep:Parameters.ComputationTStep:NbTSteps*TStep;
    subplot(4,1,i)
    plot(xis,Data.BuiltTraj(:,9))
    ylim([0 1])
    title('Condom Use')
    xlabel('time')
end
    
