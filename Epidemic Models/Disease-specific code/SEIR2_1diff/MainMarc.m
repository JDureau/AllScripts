% HPA data

DataPath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/HPA';

A = load([DataPath '/andre_estimates_31_01.txt']);

Data.Dates = {};
Data.NewCases = {};
for i = 1:size(A,1)
    Data.Dates{i} = i;
    for j = 1:7
    Data.NewCases{i}{j} = A(i,j)*10;
    end
end

InitialDate = struct();
InitialDate.Month = 6;
InitialDate.Day = 1;
InitialDate.Year = 2009;
Data.Dates = ApproxWeeklyDates(InitialDate,35);

plot(A)
legend('0-4','5-14','15-24','25-44','45-64','65+')

PopWeigths = [667600,2461800,5904100,6862500,14417400,12847800,7929300];
PopProps = PopWeigths/sum(PopWeigths)*100;
Weigthed = sum((A(:,1:7)*diag(PopWeigths.^-1)*diag(PopProps/100)*100000)');
Students = sum((A(:,1:4)*diag(PopWeigths(1:4).^-1)*diag(PopProps(1:4)/100)*100000)');
Adults = sum((A(:,5:7)*diag(PopWeigths(5:7).^-1)*diag(PopProps(5:7)/100)*100000)');


plot(A*diag(PopWeigths.^-1)*100000)
hold on
plot(Weigthed,'k','LineWidth',2)
% plot(Weigthed2,'--k','LineWidth',2)
hold off
legend('<1','1-4','5-14','15-24','25-44','45-64','65+')

for i =1:length(Data.Dates)
    disp(i)
    disp(Data.Dates{i})
end




%% Inference
cd('/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllScripts/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/SEIR'])


% 
% 
% cd('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Filtering')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\General Tools')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Parameter Estimation')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Joint Sampling')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\MIF')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Model Selection')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Optimization Approach')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\Models')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes')
% addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\SEIR')




Data.Observations = zeros(7,35);
Data.Observations(5,:) = Weigthed*10;

DiffType = 'Add';
ObsType = 'Estimated';

SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';


Name = [SavePath '/MarcData_Add_LogNorm_Tau3.mat'];

FullSEIRinference(Data,Difftype,Obstype,Name)



Thetas = [];
TransfThetas = [];
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';

inds = [1:8 10:14 16:20];
for i = 1:20%1:inds%randperm(length(inds))
    load([SavePath '/MarcCluster' num2str((i)) '.mat'])
    Thetas = [Thetas Res3.Thetas];
    TransfThetas = [TransfThetas Res3.TransfThetas];
end

Names = Res3.Parameters.Names.Estimated;
for i = 1:8
    disp(Names{i})
    ind = Res3.Parameters.(Names{i}).Index;
    disp(mean(Thetas(ind,:)))
    disp(median(Thetas(ind,:)))
    disp(quantile(Thetas(ind,:),0.025))
    disp(quantile(Thetas(ind,:),0.975))
end
    
    
    
Paths = Res3.Paths;
Parameters = Res3.Parameters;

clf
Resol = 5;
MyGrey = [237 237 237]*255^-1;
for k = 1:1
    
    %%%%%%%%%%%ResAdd
    Res = Res3;
    Parameters = Res.Parameters;
    Data = Res.Data;
    Paths = Res.Paths;
    
    dates = {};
    delta = floor(length(Res.Data.Dates)/Resol);
    Resol = floor(length(Res.Data.Dates)/delta);
   
    inds = delta:delta:delta*Resol;
    inds = [8 15 22 29];
    
    for i = 1:length(inds)
        dates{i} = [num2str(Res.Data.Dates{inds(i)}.Day) ' ' Res.Data.Dates{inds(i)}.MonthInLetters];
    end
    
    subplot(2,1,1)
    
    toplot = [5];
    for i = 1:length(toplot)
        ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.025),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.975),[172,215,255]/255)
        hold on
        ymax = 2000;
        h = area([7 13],[ymax ymax],'FaceColor',MyGrey,'EdgeColor',MyGrey);
        h = area([20 22],[ymax ymax],'FaceColor',MyGrey,'EdgeColor',MyGrey);
        set(gca,'Layer','top')
        ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.025),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.975),[172,215,255]/255)
        ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.25),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.75),[100,153,251]/255)
        plot(mean(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps))))),'k','LineWidth',2)
        plot(Data.Observations(5,:),'g','LineWidth',2)
        hold off
        xlim([1 28])
        inds = [ 7 14 21 28]
        set(gca,'XTick',inds)
%         set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
        set(gca,'XTickLabel',dates)

        if k == 1
            title('HPA estimates of total influenza incidence')
        end
    end
    
    subplot(2,1,2)
    

    try
        temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    catch
        temp = diag(Res.Thetas(Parameters.gammam1.Index,:))*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    end
    nsmooth = 10;
    temp = squeeze(exp(Paths(:,6,:)));
    ciplot(smooth(quantile(temp,0.025),nsmooth),smooth(quantile(temp,0.975),nsmooth),[172,215,255]/255)
    hold on
    h = area([Data.Instants(8-1) 252],[ymax ymax],'FaceColor',MyGrey,'EdgeColor',MyGrey);
    h = area([Data.Instants(21-1) 441],[ymax ymax],'FaceColor',MyGrey,'EdgeColor',MyGrey);
    set(gca,'Layer','top')
    ciplot(smooth(quantile(temp,0.025),nsmooth),smooth(quantile(temp,0.975),nsmooth),[172,215,255]/255)
    ciplot(smooth(quantile(temp,0.25),nsmooth),smooth(quantile(temp,0.75),nsmooth),[100,153,251]/255)
    plot(mean(temp),'k','LineWidth',2)
    t = Data.Instants-Data.NbComputingSteps(end);
    hold off
    xlim([0 568])
    TicksInds = [21*6 21*13 21*20 21*27]+1;
    set(gca,'XTick',TicksInds)
    set(gca,'XTickLabel',dates)
    if k == 1
        title('\beta_t "real-time" estimates')
    end
    ymax = 3;
    ylim([0.5 ymax])
    hold on
    yis1 = ymax/100:ymax/100:ymax*2/5;
    yis2 = ymax*3/5:ymax/100:ymax;

    hold off

    
end
    

Thetas = [];
TransfThetas = [];

inds = [1:4 6 8 10 12:20];
for i = 2:20%randperm(length(inds))
    load([SavePath '/MarcCluster' num2str((i)) '.mat'])
    Thetas = [Thetas Res3.Thetas];
    TransfThetas = [TransfThetas Res3.TransfThetas];
end
clf
plot(Thetas(1,:))



titles = {'\tau','k^{-1}','\gamma^{-1}','\beta(0)','E(0)','I(0)','R(0)','\sigma'}
clf
for i = 1:8
   subplot(3,3,i)
   plot(Thetas(i,:))
   title(titles{i})
    xlim([1 100000])
    m = mean(Thetas(i,:));
    st = std(Thetas(i,:));
    ylim([m-5*st m+5*st])
end

load([SavePath '/MarcCluster' num2str((21)) '.mat'])


Names = Res3.Parameters.Names.Estimated;
clf
for i = 1:8
   subplot(3,3,i)
   [fi,xi] = ksdensity(Thetas(i,:));
   area(xi,smooth(fi,10))
   hold on
   fis2 = zeros(size(fi));
   for j = 1:length(fis2)
       fis2(j) = normpdf(xi(j),Res3.Parameters.(Names{i}).MeanPrior,Res3.Parameters.(Names{i}).StdPrior);
   end
   if Res3.Parameters.(Names{i}).StdPrior <100
       fis2 = fis2*max(fi)/max(fis2);
%    else
%        fis2 = fis2*max(fi)/max(fis2)*(xi(end)-xi(1))*(Res3.Parameters.(Names{i}).MaxLim-Res3.Parameters.(Names{i}).MinLim);
   end
   plot(xi,fis2,'g','LineWidth',2)
   hold off
   title(titles{i})
%     xlim([1 80000])
    
end
legend('Posterior marginal density','Prior marginal density')

SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/ResultsMarc';

load([SavePath '/MarcCluster' num2str((21)) '.mat'])




days = [7 14 21 28];
for i = 1:length(days)
    Name = [SavePath 'TestsRealTime_IBM_' num2str(days(i)) 'weeks.mat'];
    Difftype = 'IBM';
    ObsType = 'Fixed';
    Data.Observations = zeros(7,days(i));
    Data.Observations(5,:) = Weigthed(1:days(i))*10;
    
    FullSEIRinference(Data,Difftype,ObsType,Name)
end
    
    
Name = [SavePath 'MarcData_IBM.mat'];

SavePath = 'S:\Results\';
load([SavePath 'ParametersSEIR.mat'])

Parameters.DiffusionType = 'Add';
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;




Difftype = 'Add';

Coeff = [5 10  20 30 40 50 60 70 80 90 150];
days = [14];
for i = 6:7%length(Coeff)
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsCoeffEst' num2str(Coeff(i)) '_Add_' num2str(days(1)) 'weeks.mat'];
    Difftype = 'Add';
    ObsType = 'CoeffStudy';
    Data.Observations = zeros(7,days(1));
    Data.Observations(5,:) = Weigthed(1:days(1))*10;
    
    Data.Coeff = Coeff(i);
    
    FullSEIRinference(Data,Difftype,ObsType,Name)
end

Difftype = 'Add';

Coeff = [5 10  20 30 40 50 70 80 ];
days = [14];
days = [14];
for i = 1:length(Coeff)
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsCoeffEst' num2str(Coeff(i)) '_Add_' num2str(days(1)) 'weeks.mat'];
    load(Name)
    PlotMarc(Res3,5)
    title(Coeff(i))
    pause()
end

Ress = {};
Coeff = [5 10  20 30 40 50 70 80 90 150];
days = [14];
days = [14];
for i = 1:length(Coeff)
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsCoeffEst' num2str(Coeff(i)) '_Add_' num2str(days(1)) 'weeks.mat'];
    load(Name)
    Res3.Coeff = Coeff(i);
    Ress{i} = Res3;
end

Diffs = [];
Means = [];
quant25 = [];
quant75 = [];
quant025 = [];
quant975 = [];
for i = 1:length(Coeff)
    Diffs(i,:) = exp(Ress{i}.Paths(:,6,150))-exp(Ress{i}.Paths(:,6,100));
    Means(i) = mean(Diffs(i,:));
    quant25(i) = quantile(Diffs(i,:),0.25);
    quant75(i) = quantile(Diffs(i,:),0.75);
    quant025(i) = quantile(Diffs(i,:),0.025);
    quant975(i) = quantile(Diffs(i,:),0.975);
end

clf
plot(Coeff,Means,'k','LineWidth',2)
hold on
jbfill(Coeff,quant975,quant025,[172,215,255]/255);
jbfill(Coeff,quant75,quant25,[100,153,251]/255);
jbfill(Coeff,Means,Means*1.05,0*[255,255,255]/255)
jbfill(Coeff,zeros(size(Coeff)),zeros(size(Coeff)),0*[255,255,255]/255)
hold off
xlim([Coeff(1) Coeff(end)])
xlabel('Multiplicative correction coefficient')
ylabel('\beta(July 13^{th}) - \beta(August 1^{st})')
legend('Mean estimate','95% CI','50% CI')


Difftype = 'Add';


days = [14];
Coeff = [5 10  40 90 150];
for i = 4:length(Coeff)
    Name = [SavePath 'TestsCoeffEst' num2str(Coeff(i)) '_Add_' num2str(days(1)) 'weeks.mat'];
    Difftype = 'Add';
    ObsType = 'CoeffEst';
    Data.Observations = zeros(7,days(1));
    Data.Observations(5,:) = Weigthed(1:days(1))*10;
    
    load(Name)
    PlotMarc(Res3,5)
end










Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

    

SEIRModel = struct();
SEIRModel.UKF_projection = @SEIR_UKF_projection;
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

temp = zeros(1,7);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end




Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters.InitialCovFact.Estimated =  1;

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


test = 0;
NbIts = 0;
while not(test)
    Parameters = SampleParameters(Parameters);
    Parameters.SigmaRW.Value = rand(1,1)*2;
    Parameters.InitialCovFact.Value = rand(1,1)*0.5;
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    try
        Temp = EstimationEKFGen(Data, SEIRModel, Parameters);
        Temp.LogLik
        if Temp.LogLik>-300
            test = 1;
        end
    end
    NbIts = NbIts + 1;
    if NbIts>10000
        test = 1;
    end
end

    

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',15000,'TolX',1e-8,'TolFun',1e-7));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  



Parameters.InitialCovFact.Value = 0.1;
Parameters = UpdateParsNoTransfToTransf(Parameters);

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;


Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);




Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',15000,'TolX',1e-8,'TolFun',1e-7));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  



Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).SamplStd = 0.0001*Parameters.(Names{i}).Value;
end
ResKal = KalmanNumericDerivativesWithPrior(Data,SEIRModel,Parameters);

Test = mean(eig(-ResKal.Hess)>0)==1;
disp(Test)
Cov = (-ResKal.Hess)^-1;

%% SMC Optimization

Parameters.ComputationTStep = 0.1;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = A(:,7);
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',100,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


% defining NbParts 
% nbparts = [200 400 600 800 1000 1200 1400 1600];
% NbTests = 100;
% LogLikrecords = [];
% for i = 1:length(nbparts)
%     Parameters.NbParticules = nbparts(i);
%     for j = 1:NbTests
%         disp([i j])
%         Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
%         LogLikrecords(i,j) = Temp.LogLik;
%     end
% end
% clf
% hold on
% for i = 1:length(nbparts)
%     plot(nbparts(i),std(LogLikrecords(i,:)),'o')
% end
% hold off
% 
% Parameters.NbParticules = 1000;
% 
% % defining TStep
% DataTests = Data;
% 
% res = [1:10].^-1;
% logliks = [];
% paths = [];
% for i = 1:length(res)
%     disp(i)
%     for j = 1:10
%         Parameters.ComputationTStep = res(i);
%         DataTests.Observations = zeros(6,35);
%         DataTests.Observations(5,:) = A(:,7);
%         DataTests.Instants = [1:35]*7/Parameters.ComputationTStep;
%         DataTests.ObservedVariables = 5*ones(1,length(DataTests.Instants));
%         DataTests.NbComputingSteps = [0 diff(DataTests.Instants)];
%         Temp = EstimationSMCsmoothGen(DataTests, SEIRModel, Parameters);
%         logliks(i,j) = Temp.LogLik;
%         paths(i,j,:) = exp(Temp.PosteriorMeansRecord(6,:));
%     end
% end
% clf
% hold on
% for i = 1:length(res)
%     plot(res(i),mean(logliks(i,:)),'.')
%     plot(res(i),mean(logliks(i,:))+std(logliks(i,:)),'.r')
%     plot(res(i),mean(logliks(i,:))-std(logliks(i,:)),'.r')
% end
% hold off
% 
% cols = rand(length(res),3);
% clf
% hold on
% for i = 1:length(res)
%     plot(mean(squeeze(paths(i,:,:))),'col',cols(i,:),'LineWidth',2)
%     plot(quantile(squeeze(paths(i,:)),0.95),'col',cols(i,:))
%     plot(quantile(squeeze(paths(i,:)),0.05),'col',cols(i,:))
% end
% legend
% hold off


Parameters.ComputationTStep = 1/3;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = Weigthed*10;
Data.Instants = [0:34]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.ModelType = 'SMC';
dim = 7;
Cov = 2.38^2/dim*(-ResKal.Hess)^-1;
Parameters.G = -ResKal.Hess;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 1;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);

TempRes = Res;

dim = length(Parameters.Names.Estimated);
Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);

SavePath = 'S:\Results\';
save([SavePath 'Marc_ForPresEst_NoPaths.mat'],'Res2')


TempRes = Res2;

Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.6;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);
Res3.Description = 'This one has unif prior on days for latent and inf periods. It also uses the actual averaged data and not only >65 as FirstEst';
SavePath = 'S:\Results\';
save([SavePath 'MarcData_Add.mat'],'Res3')

Res = load([SavePath 'Marc_ForPresEst_WithPaths.mat'])

Res3 = Res.Res3;
tmp = mean(Res3.Thetas');
Parameters = Res3.Parameters;
Names = Parameters.Names.Estimated;
for i = 1:length(tmp)
    Parameters.(Names{i}).Value = tmp(Parameters.(Names{i}).Index);
end
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters);
TempThetaMean = EstimationSMCsmoothGen(Res3.Data, SEIRModel, Parameters);

p_d = mean(-2*Res3.LogLiks) - (-2)*TempThetaMean.LogLik;
DIC = (-2)*TempThetaMean.LogLik + 2*p_d;

Res3.p_d = p_d;
Res3.DIC = DIC;

save([SavePath 'Marc_ForPresEst_WithPaths.mat'],'Res3')

%% IBM


Parameters = struct();

Parameters.Problem = 'MarcFlu';
Parameters.NbVariables = 7;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'IBM';
Parameters.ObservationLength = 7*35;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation = 100000;

Data.Observations = zeros(7,35);
Data.Observations(5,:) = Weigthed*10;
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];



Parameters.km1.Value = 1.6;
Parameters.km1.Value = 1.6;
Parameters.km1.Min = -10^14;
Parameters.km1.Max = 10^14;
Parameters.km1.MinLim = 1.55;
Parameters.km1.MaxLim = 1.63;
Parameters.km1.Estimated = 1;
Parameters.km1.TransfType = 'Logit';
Parameters.km1.PlotInv = 1;
Parameters.km1.Sample = 1;
Parameters.gammam1.Value = 1;
Parameters.gammam1.Min = -10^14;
Parameters.gammam1.Max = 10^14;
Parameters.gammam1.MinLim = 0.93;
Parameters.gammam1.MaxLim = 1.23;
Parameters.gammam1.Estimated = 1;
Parameters.gammam1.TransfType = 'Logit';
Parameters.gammam1.PlotInv = 1;
Parameters.gammam1.Sample = 1;
Parameters.betainit.Value = 0.25;
Parameters.betainit.Min = -10^14;
Parameters.betainit.Max = 10^14;
Parameters.betainit.Estimated = 1;
Parameters.betainit.TransfType = 'Log';
Parameters.betainit.Init = 1;
Parameters.betainit.Sample = 1;
Parameters.betaderinit.Value = 0;
Parameters.betaderinit.Min = -10^14;
Parameters.betaderinit.Max = 10^14;
Parameters.betaderinit.MinLim = -100;
Parameters.betaderinit.MaxLim =  100;
Parameters.betaderinit.Estimated = 0;
Parameters.betaderinit.TransfType = 'Logit';
Parameters.betaderinit.Init = 1;
Parameters.betaderinit.Sample = 0;
Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MinLim = 0;
Parameters.EInitProp.MaxLim = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.TransfType = 'Logit';
Parameters.EInitProp.Init = 1;
Parameters.EInitProp.Sample = 0;
Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MinLim = 0;
Parameters.IInitProp.MaxLim = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.IInitProp.TransfType = 'Logit';
Parameters.IInitProp.Init = 1;
Parameters.IInitProp.Sample = 0;
Parameters.RInitProp.Value = 0.2;
Parameters.RInitProp.Min = 0;
Parameters.RInitProp.Max = 0.30;
Parameters.RInitProp.MinLim = 0;
Parameters.RInitProp.MaxLim = 0.6;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.TransfType = 'Logit';
Parameters.RInitProp.Init = 1;
Parameters.RInitProp.Sample = 0;
Parameters.SigmaRW.Value = exp(-0.6);
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters.SigmaRW.Sample = 0;
Parameters.InitialCovFact.Value = 0.20;
Parameters.InitialCovFact.Min = -10^14;
Parameters.InitialCovFact.Max =  10^14;
Parameters.InitialCovFact.Estimated =  0;
Parameters.InitialCovFact.TransfType = 'Log';
Parameters.InitialCovFact.Init = 1;
Parameters.InitialCovFact.Sample = 0;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

% Parameters.MIFNbIterations = 500;
% Parameters.MIFNbParticules = 3000;
% Parameters.MIFCoolingParameters = 0.99;
% Parameters.MIFb = 1;
% Parameters.MIFSigmas = [];
% Parameters.NbParticules = 1000;
% Parameters.NoPaths = 0;
% Parameters.PathsToKeep = [1:6]';
% Parameters.MCMCType = 'gtgt';
% Names = Parameters.Names.Estimated;
% for i = 1:length(Parameters.Names.Estimated)
%     Parameters.MIFSigmas(i) = 0.1*Parameters.(Names{i}).TransfValue;
% end
% Result = MIFSEIR(Data,Parameters,SEIRModel)
% plot(Result.LogLiks)
% 
% TrueParameters = Res.Res3.Parameters;
% k = Parameters.NbParsEstimated;
% for i = 1:k
%     subplot(k,1,i)
%     plot(Result.ThetasRecord(i,:))
%     hold on
%     xs = 1:length(Result.ThetasRecord(i,:));
%     plot(xs,mean(Res.Res3.TransfThetas(i,4000:5000)),'r')
%     hold off
% end
% 
% Temp = EstimationSMCsmoothGen(Data, SEIRModel, TrueParameters);
% PlotMarc(Temp,8)
% 
%     
    

SEIRModel = struct();
SEIRModel.UKF_projection = @SEIR_UKF_projection;
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

% EKF Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];



temp = zeros(1,7);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end



Res = load([SavePath 'Marc_ForPresEst_WithPaths.mat']);
Res = Res.Res3;
ParametersForCov = Res.Parameters;
tmp = mean(Res.Thetas');
Names = ParametersForCov.Names.Estimated;
for i = 1:length(tmp)
    ParametersForCov.(Names{i}).Value = tmp(ParametersForCov.(Names{i}).Index);
end
ParametersForCov = UpdateParsNoTransfToTransf(ParametersForCov);


Names = ParametersForCov.Names.Estimated;
for i = 1:length(ParametersForCov.Names.Estimated)
    Parameters.(Names{i}).Value = ParametersForCov.(Names{i}).Value;
end
Parameters = UpdateParsNoTransfToTransf(Parameters);
    


Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.betaderinit.Estimated = 1;
Parameters.InitialCovFact.Estimated = 1;
Parameters.InitialCovFact.Value = 0.3;
Parameters.gammam1.Estimated = 1;
Parameters.km1.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);


Test = 0;
NbIts = 0;
while not(Test)
    Parameters.SigmaRW.Value = rand(1,1)*2;
    Parameters.InitialCovFact.Value = rand(1,1)*0.5;
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    try
        Temp = EstimationEKFGen(Data, SEIRModel, Parameters);
        if (Temp.LogLik>-1000)
            Test = 1;
        end
    end
    NbIts = NbIts + 1;
    if NbIts>20000
        disp('Can''t initialize IBM')
        die
    end
end
    

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',15000,'TolX',1e-8,'TolFun',1e-7));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end

Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters.InitialCovFact.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);




Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',15000,'TolX',1e-8,'TolFun',1e-7));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).SamplStd = 0.01*Parameters.(Names{i}).Value;
end
ResKal = KalmanNumericDerivativesWithPrior(Data,SEIRModel,Parameters);

Test = mean(eig(-ResKal.Hess)>0)==1;
disp(Test)
Cov = (-ResKal.Hess)^-1;

%% SMC Optimization



Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',100,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


% defining NbParts 
% nbparts = [200 400 600 800 1000 1200 1400 1600];
% NbTests = 100;
% LogLikrecords = [];
% for i = 1:length(nbparts)
%     Parameters.NbParticules = nbparts(i);
%     for j = 1:NbTests
%         disp([i j])
%         Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
%         LogLikrecords(i,j) = Temp.LogLik;
%     end
% end
% clf
% hold on
% for i = 1:length(nbparts)
%     plot(nbparts(i),std(LogLikrecords(i,:)),'o')
% end
% hold off
% 
% Parameters.NbParticules = 1000;
% 
% % defining TStep
% DataTests = Data;
% 
% res = [1:10].^-1;
% logliks = [];
% paths = [];
% for i = 1:length(res)
%     disp(i)
%     for j = 1:10
%         Parameters.ComputationTStep = res(i);
%         DataTests.Observations = zeros(6,35);
%         DataTests.Observations(5,:) = A(:,7);
%         DataTests.Instants = [1:35]*7/Parameters.ComputationTStep;
%         DataTests.ObservedVariables = 5*ones(1,length(DataTests.Instants));
%         DataTests.NbComputingSteps = [0 diff(DataTests.Instants)];
%         Temp = EstimationSMCsmoothGen(DataTests, SEIRModel, Parameters);
%         logliks(i,j) = Temp.LogLik;
%         paths(i,j,:) = exp(Temp.PosteriorMeansRecord(6,:));
%     end
% end
% clf
% hold on
% for i = 1:length(res)
%     plot(res(i),mean(logliks(i,:)),'.')
%     plot(res(i),mean(logliks(i,:))+std(logliks(i,:)),'.r')
%     plot(res(i),mean(logliks(i,:))-std(logliks(i,:)),'.r')
% end
% hold off
% 
% cols = rand(length(res),3);
% clf
% hold on
% for i = 1:length(res)
%     plot(mean(squeeze(paths(i,:,:))),'col',cols(i,:),'LineWidth',2)
%     plot(quantile(squeeze(paths(i,:)),0.95),'col',cols(i,:))
%     plot(quantile(squeeze(paths(i,:)),0.05),'col',cols(i,:))
% end
% legend
% hold off
% 

Parameters.ComputationTStep = 1/3;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = Weigthed*10;
Data.Instants = [0:34]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.ModelType = 'SMC';
dim = 7;
Cov = 2.38^2/dim*(-ResKal.Hess)^-1;
Parameters.G = -ResKal.Hess;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 1;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,100);

TempRes = Res;

dim = length(Parameters.Names.Estimated);
Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);

SavePath = 'S:\Results\';
save([SavePath 'Marc_ForPresEst_NoPaths_IBM.mat'],'Res2')


TempRes = Res2;

Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.6;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);
Res3.Description = 'This one has unif prior on days for latent and inf periods. It also uses the actual averaged data and not only >65 as FirstEst';
SavePath = 'S:\Results\';
save([SavePath 'Marc_ForPresEst_WithPaths.mat'],'Res3')

Res = load([SavePath 'Marc_ForPresEst_WithPaths_IBM.mat'])






%% VolVol


Parameters = struct();

Parameters.Problem = 'MarcFlu';
Parameters.NbVariables = 7;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'SVO';
Parameters.ObservationLength = 7*35;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation = 100000;

Data.Observations = zeros(7,35);
Data.Observations(5,:) = Weigthed*10;
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];



Parameters.km1.Value = 1.6;
Parameters.km1.Min = -10^14;
Parameters.km1.Max = 10^14;
Parameters.km1.MinLim = 1.55;
Parameters.km1.MaxLim = 1.63;
Parameters.km1.Estimated = 1;
Parameters.km1.TransfType = 'Logit';
Parameters.km1.PlotInv = 1;
Parameters.gammam1.Value = 1;
Parameters.gammam1.Min = -10^14;
Parameters.gammam1.Max = 10^14;
Parameters.gammam1.MinLim = 0.93;
Parameters.gammam1.MaxLim = 1.23;
Parameters.gammam1.Estimated = 1;
Parameters.gammam1.TransfType = 'Logit';
Parameters.gammam1.PlotInv = 1;
Parameters.betainit.Value = 0.25;
Parameters.betainit.Min = -10^14;
Parameters.betainit.Max = 10^14;
Parameters.betainit.Estimated = 1;
Parameters.betainit.TransfType = 'Log';
Parameters.betainit.Init = 1;
Parameters.betaderinit.Value = 0;
Parameters.betaderinit.Min = -10^14;
Parameters.betaderinit.Max = 10^14;
Parameters.betaderinit.Estimated = 0;
Parameters.betaderinit.TransfType = 'Log';
Parameters.VolInit.Init = 1;
Parameters.VolInit.Value = 0;
Parameters.VolInit.Min = -10^14;
Parameters.VolInit.Max = 10^14;
Parameters.VolInit.Estimated = 1;
Parameters.VolInit.TransfType = 'Log';
Parameters.VolInit.Init = 1;
Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.EInitProp.MinLim = 0;
Parameters.EInitProp.MaxLim = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.EInitProp.TransfType = 'Logit';
Parameters.EInitProp.Init = 1;
Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Min = 0.002*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Max = 10*max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.MinLim = 0;
Parameters.IInitProp.MaxLim = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.IInitProp.TransfType = 'Logit';
Parameters.IInitProp.Init = 1;
Parameters.RInitProp.Value = 0.2;
Parameters.RInitProp.Min = 0;
Parameters.RInitProp.Max = 0.30;
Parameters.RInitProp.MinLim = 0;
Parameters.RInitProp.MaxLim = 0.6;
Parameters.RInitProp.Estimated = 1;
Parameters.RInitProp.TransfType = 'Logit';
Parameters.RInitProp.Init = 1;
Parameters.SigmaRW.Value = exp(-0.6);
Parameters.SigmaRW.Min = -10^14;
Parameters.SigmaRW.Max = 10^14;
Parameters.SigmaRW.Estimated = 1;
Parameters.SigmaRW.TransfType = 'Log';
Parameters.InitialCovFact.Value = 0.20;
Parameters.InitialCovFact.Min = -10^14;
Parameters.InitialCovFact.Max =  10^14;
Parameters.InitialCovFact.Estimated =  2;
Parameters.InitialCovFact.TransfType = 'Log';
Parameters.InitialCovFact.Init = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

% Parameters.MIFNbIterations = 500;
% Parameters.MIFNbParticules = 3000;
% Parameters.MIFCoolingParameters = 0.99;
% Parameters.MIFb = 1;
% Parameters.MIFSigmas = [];
% Parameters.NbParticules = 1000;
% Parameters.NoPaths = 0;
% Parameters.PathsToKeep = [1:6]';
% Parameters.MCMCType = 'gtgt';
% Names = Parameters.Names.Estimated;
% for i = 1:length(Parameters.Names.Estimated)
%     Parameters.MIFSigmas(i) = 0.1*Parameters.(Names{i}).TransfValue;
% end
% Result = MIFSEIR(Data,Parameters,SEIRModel)
% plot(Result.LogLiks)
% 
% TrueParameters = Res.Res3.Parameters;
% k = Parameters.NbParsEstimated;
% for i = 1:k
%     subplot(k,1,i)
%     plot(Result.ThetasRecord(i,:))
%     hold on
%     xs = 1:length(Result.ThetasRecord(i,:));
%     plot(xs,mean(Res.Res3.TransfThetas(i,4000:5000)),'r')
%     hold off
% end
% 
% Temp = EstimationSMCsmoothGen(Data, SEIRModel, TrueParameters);
% PlotMarc(Temp,8)
% 
%     
    

SEIRModel = struct();
SEIRModel.UKF_projection = @SEIR_UKF_projection;
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),Data.Observations(5,IndTime),Data.Observations(5,IndTime)*Parameters.SigmaObs)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);

% EKF Optimization
Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];



temp = zeros(1,7);
temp(1,5) = 1;
SEIRModel.ObservationJacobian = {};
SEIRModel.ObservationMeasurementNoise = {};
for i = 1:length(Data.Instants)
    SEIRModel.ObservationJacobian{i} = temp;
    SEIRModel.ObservationMeasurementNoise{i} = (Parameters.SigmaObs*Data.Observations(5,i))^2;
end



Res = load([SavePath 'Marc_ForPresEst_WithPaths.mat']);
Res = Res.Res3;
Parameters = Res.Parameters;

Parameters.DiffusionType = 'SVO';

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.SigmaRW.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.VolInit.Init = 0.2;
Parameters.VolInit.Value = 0.5;
Parameters.VolInit.Min = -10^14;
Parameters.VolInit.Max = 10^14;
Parameters.VolInit.Estimated = 1;
Parameters.VolInit.TransfType = 'Log';
Parameters.VolInit.Init = 1;
Parameters.InitialCovFact.Estimated = 1;
Parameters.InitialCovFact.Value = 0.3;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);




Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',15000,'TolX',1e-8,'TolFun',1e-7));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end

Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters.InitialCovFact.Estimated = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);




Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) KalmanToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',15000,'TolX',1e-8,'TolFun',1e-7));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).SamplStd = 0.01*Parameters.(Names{i}).Value;
end
ResKal = KalmanNumericDerivativesWithPrior(Data,SEIRModel,Parameters);

Test = mean(eig(-ResKal.Hess)>0)==1;
disp(Test)
Cov = (-ResKal.Hess)^-1;

%% SMC Optimization

Parameters.ComputationTStep = 0.1;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = A(:,7);
Data.Instants = [1:35]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.NoPaths = 1;
Parameters.NbParticules = 1000;
Parameters.MCMCType = 'frfr';
Names = Parameters.Names.Estimated;
Initialization = [];

Names = Parameters.Names.Estimated;
for j = 1:length(Names)
    Parameters.(Names{j}).Estimated = 0;
end
Parameters.km1.Estimated = 1;
Parameters.gammam1.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.EInitProp.Estimated = 1;
Parameters.IInitProp.Estimated = 1;
Parameters.RInitProp.Estimated = 1;
Parameters.SigmaRW.Estimated = 1;
Parameters.betainit.Estimated = 1;
Parameters.VolInit.Init = 1;
Parameters.VolInit.Value = 0.5;
Parameters.VolInit.Min = -10^14;
Parameters.VolInit.Max = 10^14;
Parameters.VolInit.Estimated = 1;
Parameters.VolInit.TransfType = 'Log';
Parameters.VolInit.Init = 1;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Initialization(i) = Parameters.(Names{i}).TransfValue ;
%         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
end
SEIRModel.InitializeParameters = @SEIRInitialize;
[x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',100,'TolX',1e-8,'TolFun',1e-8));
Names = Parameters.Names.Estimated;
for i = 1:length(Names)
    Parameters.(Names{i}).TransfValue = (x(i));
end
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)  


% defining NbParts 
nbparts = [200 400 600 800 1000 1200 1400 1600];
NbTests = 100;
LogLikrecords = [];
for i = 1:length(nbparts)
    Parameters.NbParticules = nbparts(i);
    for j = 1:NbTests
        disp([i j])
        Temp = EstimationSMCsmoothGen(Data, SEIRModel, Parameters);
        LogLikrecords(i,j) = Temp.LogLik;
    end
end
clf
hold on
for i = 1:length(nbparts)
    plot(nbparts(i),std(LogLikrecords(i,:)),'o')
end
hold off

Parameters.NbParticules = 1000;

% defining TStep
DataTests = Data;

res = [1:10].^-1;
logliks = [];
paths = [];
for i = 1:length(res)
    disp(i)
    for j = 1:10
        Parameters.ComputationTStep = res(i);
        DataTests.Observations = zeros(6,35);
        DataTests.Observations(5,:) = A(:,7);
        DataTests.Instants = [1:35]*7/Parameters.ComputationTStep;
        DataTests.ObservedVariables = 5*ones(1,length(DataTests.Instants));
        DataTests.NbComputingSteps = [0 diff(DataTests.Instants)];
        Temp = EstimationSMCsmoothGen(DataTests, SEIRModel, Parameters);
        logliks(i,j) = Temp.LogLik;
        paths(i,j,:) = exp(Temp.PosteriorMeansRecord(6,:));
    end
end
clf
hold on
for i = 1:length(res)
    plot(res(i),mean(logliks(i,:)),'.')
    plot(res(i),mean(logliks(i,:))+std(logliks(i,:)),'.r')
    plot(res(i),mean(logliks(i,:))-std(logliks(i,:)),'.r')
end
hold off

cols = rand(length(res),3);
clf
hold on
for i = 1:length(res)
    plot(mean(squeeze(paths(i,:,:))),'col',cols(i,:),'LineWidth',2)
    plot(quantile(squeeze(paths(i,:)),0.95),'col',cols(i,:))
    plot(quantile(squeeze(paths(i,:)),0.05),'col',cols(i,:))
end
legend
hold off


Parameters.ComputationTStep = 1/3;

Data.Observations = zeros(6,35);
Data.Observations(5,:) = Weigthed*10;
Data.Instants = [0:34]*7/Parameters.ComputationTStep;
Data.ObservedVariables = 5*ones(1,length(Data.Instants));
Data.NbComputingSteps = [0 diff(Data.Instants)];


Parameters.ModelType = 'SMC';
dim = 7;
Cov = 2.38^2/dim*(-ResKal.Hess)^-1;
Parameters.G = -ResKal.Hess;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 1;
TempPar = ProposeInitialParameter(Data, SEIRModel, Parameters);
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);

TempRes = Res;

dim = length(Parameters.Names.Estimated);
Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.NoPaths = 1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.8;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Res2 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,5000);

SavePath = 'S:\Results\';
save([SavePath 'Marc_ForPresEst_NoPaths.mat'],'Res2')


TempRes = Res2;

Cov = 2.38^2/dim*cov(TempRes.TransfThetas');
Parameters.G = Cov^-1;
Parameters.MCMCType = 'Rand';
Parameters.GMeth = 'cst given';
Parameters.ComputeRWsamplCov = 0;
Parameters.aim = 0.23;
Parameters.Epsil = 0.6;
TempPar = TempRes.TempPar;
% [Parameters, TempPar] = CalibrateMethod( Data, SEIRModel, Parameters, TempPar);
Parameters.NoPaths = 0;
Parameters.PathsToKeep = [1:7]';
Res3 = RunEstimationMethod(Data, SEIRModel,Parameters,TempPar,10000);
Res3.Description = 'This one has unif prior on days for latent and inf periods. It also uses the actual averaged data and not only >65 as FirstEst';
SavePath = 'S:\Results\';
save([SavePath 'Marc_ForPresEst_WithPaths.mat'],'Res3')

Res = load([SavePath 'Marc_ForPresEst_WithPaths_IBM.mat'])
















