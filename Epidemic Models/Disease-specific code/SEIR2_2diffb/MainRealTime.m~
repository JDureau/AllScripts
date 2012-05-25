%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           DOES IT WORK IN REAL TIME?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% HPA data

DataPath = 'H:\My Documents\PhD Work\Data\HPA';

A = load([DataPath '\andre_estimates_31_01.txt']);

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
cd('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Filtering')
addpath('H:\My Documents\PhD Work\Matlab Scripts\General Tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Joint Sampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\MIF')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Model Selection')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Optimization Approach')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\Models')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\SEIR')

Parameters = struct();

Parameters.Problem = 'MarcFlu';
Parameters.NbVariables = 7;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'Add';
Parameters.ObservationLength = 7*35;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation = 100000;

tmp = 1:4;
lengths = ceil(35/5*tmp);

for i = 1:length(lengths)
    Data.Observations = zeros(7,lengths(i));
    Data.Observations(5,:) = Weigthed(1:lengths(i))*10;
    Data.Instants = [1:lengths(i)]*7/Parameters.ComputationTStep;
    Data.ObservedVariables = 5*ones(1,length(Data.Instants));
    Data.NbComputingSteps = [0 diff(Data.Instants)];

    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsRealTime_Add_' num2str(lengths(i)) 'weeks.mat'];
    Difftype = 'Add';
    ObsType = 'Fixed';
    FullSEIRinference(Data,Difftype,ObsType,Name)
end




SavePath = 'S:\Results\';
tmp = load([SavePath 'ParametersSEIR.mat']);
Parameters = tmp.Parameters;

Parameters.EInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters.IInitProp.Value = max(1,Data.Observations(5,1))/Parameters.TotalPopulation;
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters)

SEIRModel = struct();
SEIRModel.EKF_projection = @SEIR_EKF_projection;
SEIRModel.InitializeParameters = @SEIRInitialize;
SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
SEIRModel.SMC_projection = @SEIR_SMC_projection;
Parameters = SEIRModel.InitializeParameters(Parameters);


for i = 1:length(lengths)
    Data.Observations = zeros(7,lengths(i));
    Data.Observations(5,:) = Weigthed(1:lengths(i))*10;
    Data.Instants = [1:lengths(i)]*7/Parameters.ComputationTStep;
    Data.ObservedVariables = 5*ones(1,length(Data.Instants));
    Data.NbComputingSteps = [0 diff(Data.Instants)];

    
    
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsRealTime_Add_' num2str(lengths(i)) 'weeks.mat'];
    Res = load(Name);
    
    Parameters = Res.Res3.Parameters;
    
    Parameters.MIFCov = cov(Res.Res3.TransfThetas');
    
    
    Parameters.ComputationTStep = 1/3;

    Data.Instants = [0:size(Data.Observations,2)-1]*7/Parameters.ComputationTStep;
    Data.ObservedVariables = 5*ones(1,length(Data.Instants));
    Data.NbComputingSteps = [0 diff(Data.Instants)];

    
    SEIRModel = struct();
    SEIRModel.EKF_projection = @SEIR_EKF_projection;
    SEIRModel.InitializeParameters = @SEIRInitialize;
    SEIRModel.LikFunction = 'normpdf(Variables(:,5),coeff*Data.Observations(5,IndTime),coeff*Data.Observations(5,IndTime)*Parameters.SigmaObs.Value)';
    SEIRModel.SMC_projection = @SEIR_SMC_projection;
    Parameters = SEIRModel.InitializeParameters(Parameters);

    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        Parameters.(Names{j}).Value = mean(Res.Res3.Thetas(j,:));
    end
    Parameters = DefineEstimatedParametersIndexes(Parameters);
    Parameters = DefineTransfFunctions(Parameters);
    Parameters = DefinePriors(Parameters);
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    
    Parameters.MIFNbIterations = 100;
    Parameters.MIFNbParticules = 10000;   
    Parameters.MIFCoolingParameters = 0.95;
    Parameters.MIFb= 1;
    
%     Result = MIFSEIR(Data,Parameters,SEIRModel);
%     ParamsMIF = Result.Parameters;
        
    


    Initialization = [];
    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        Initialization(j) = Parameters.(Names{j}).TransfValue ;
    %         Initialization(i) =  tempPars.(Names{i}).TransfValue ;
    end
    SEIRModel.InitializeParameters = @SEIRInitialize;
    Parameters.Correction = 0;
    [x,fval,exitflag,output] = fminsearch(@(x) SMCToOptimizeWithPrior(x,Data,SEIRModel,Parameters),Initialization,optimset('MaxIter',200,'TolX',1e-8,'TolFun',1e-8));
    Names = Parameters.Names.Estimated;
    for j = 1:length(Names)
        Parameters.(Names{j}).TransfValue = (x(j));
    end
    Parameters = UpdateParsTransfToNoTransf(Parameters);
    TellParsValues(Parameters)  
    ParamsMIF = Parameters;

%     
%     Names = Parameters.Names.Estimated;
%     xs = [];
%     LogLiks = [];
%     Cov = Parameters.MIFCov;
%     for j = 1:length(Names)
%         TempPars = ParamsMIF;
%         Result = EstimationSMCsmoothGen(Data, SEIRModel, TempPars);
%         xs(j,1) = TempPars.(Names{j}).Value;
%         LogLiks(j,1) = Result.LogLik+ Result.LogPrior;
%        
%         for k = 1:10
%             TempPars = ParamsMIF;
%             TempPars.(Names{j}).TransfValue = TempPars.(Names{j}).TransfValue + k/10*sqrt(Cov(j,j));
%             TempPars = UpdateParsTransfToNoTransf(TempPars);
%             Result = EstimationSMCsmoothGen(Data, SEIRModel, TempPars);
%             xs(j,(k-1)*2+2) = TempPars.(Names{j}).Value;
%             LogLiks(j,(k-1)*2+2) = Result.LogLik + Result.LogPrior;
%             TempPars = ParamsMIF;
%             TempPars.(Names{j}).TransfValue = TempPars.(Names{j}).TransfValue - k/10*sqrt(Cov(j,j));
%             TempPars = UpdateParsTransfToNoTransf(TempPars);
%             Result = EstimationSMCsmoothGen(Data, SEIRModel, TempPars);
%             xs(j,(k-1)*2+3) = TempPars.(Names{j}).Value;
%             LogLiks(j,(k-1)*2+3) = Result.LogLik + Result.LogPrior;
%         end    
%         plot(xs(j,:),LogLiks(j,:),'.')
%         pause()
%     end
%     

    
    SMCPaths = [];
    NbIts = 500;
    ParamesMIF.NoPaths = 0;
    ParamsMIF.PathsToKeep = [1:7]';
    ParamsMIF.NbParticules = 1000;
    for j = 1:NbIts
       j
        Temp = EstimationSMCsmoothGen(Data, SEIRModel, ParamsMIF);
        ind = ceil(rand(1,1)*ParamsMIF.NbParticules);
        SMCPaths(end+1,:,:) = Temp.CompletePaths(ind,:,:);
    end
    
    Res = struct();
    Res.SMCPaths = SMCPaths;
    Res.ParamsMIF = ParamsMIF; 
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsRealTime_MIF_' num2str(lengths(i)) 'weeks_SecondTrial.mat'];
    Difftype = 'Add';
    ObsType = 'Fixed';
    save(Name,'Res')
%     FullSEIRinference(Data,Difftype,ObsType,Name)
end




tmp = 1:4;
lengths = ceil(35/5*tmp);
RessMIF = {};
RessAdd = {};
for i = 1:length(lengths)
    
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsRealTime_MIF_' num2str(lengths(i)) 'weeks_SecondTrial.mat'];
    Res = load(Name);
    RessMIF{i} = Res.Res;
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsRealTime_Add_' num2str(lengths(i)) 'weeks.mat'];
    Res = load(Name);
    RessAdd{i} = Res.Res3;
%     figure
%     PlotMarc(Res,8)
%     pause()
end
 
for i = 1:length(Res.Res3.Parameters.Na
    mes.Estimated)
    MIFS = [];
    for j = 1:length(RessAdd)
        MIF(j) = ResMIF{j}.ParamsMIF.();
        
    end
end







for i = 1:length(lengths)
    
    
    Res = Ress{i} ;
    Paths = Res.SMCPaths;
    subplot(4,1,i)
    ciplot(quantile(exp(squeeze(Paths(:,6,:))),0.025),quantile(exp(squeeze(Paths(:,6,:))),0.975),[172,215,255]/255)
    hold on
    ciplot(quantile(exp(squeeze(Paths(:,6,:))),0.25),quantile(exp(squeeze(Paths(:,6,:))),0.75),[100,153,251]/255)
    plot(mean(exp(squeeze(Paths(:,6,:)))),'k','LineWidth',2)
    hold off
    title(Res.ParamsMIF.RInitProp.Value)
%     figure
%     PlotMarc(Res,8)
%     pause()
end



for i = 1:length(lengths)
    
    
    Res = Ress{i} ;
    Paths = Res.SMCPaths;
    subplot(4,1,i)
    tmp = exp(squeeze(Paths(:,6,:))).*(squeeze(Paths(:,1,:)))*Res.ParamsMIF.gammam1.Value;
    ciplot(quantile(tmp,0.025),quantile(tmp,0.975),[172,215,255]/255)
    hold on
    ciplot(quantile(tmp,0.25),quantile(tmp,0.75),[100,153,251]/255)
    plot(mean(tmp),'k','LineWidth',2)
    hold off
    title(Res.ParamsMIF.gammam1.Value)
    xlim([0 600])
%     figure
%     PlotMarc(Res,8)
%     pause()
end


clf
Resol = 5;
grey = [237 237 237]*255^-1;
for k = 1:4
    
    %%%%%%%%%%%ResAdd
    Res = RessAdd{k};
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
    
    subplot(4,2,2*(k-1)+1)
    
    toplot = [5];
    for i = 1:length(toplot)
        ciplot(quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.025),quantile(squeeze(Paths(:,toplot(i),max(1,cumsum(Data.NbComputingSteps)))),0.975),[172,215,255]/255)
        hold on
        ymax = 2000;
        h = area([7 13],[ymax ymax],'FaceColor',rey,'EdgeColor',grey);
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
    
        subplot(4,2,2*(k-1)+2)
    
    try
        temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    catch
        temp = diag(Res.Thetas(Parameters.gammam1.Index,:))*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
    end
    temp = squeeze(exp(Paths(:,6,:)));
    ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
    hold on
    h = area([Data.Instants(8-1) 252],[ymax ymax],'FaceColor',MyGrey,'EdgeColor',MyGrey);
    set(gca,'Layer','top')
    ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
    ciplot(quantile(temp,0.25),quantile(temp,0.75),[100,153,251]/255)
    plot(mean(temp),'k','LineWidth',2)
    t = Data.Instants-Data.NbComputingSteps(end);
%     plot(t,1*ones(size(Data.Instants)),'--k','LineWidth',2)
%     try
%         plot(Res.RtPath,'g','LineWidth',2) 
%     end
    hold off
    xlim([0 568])
    TicksInds = [21*6 21*13 21*20 21*27]+1;
    set(gca,'XTick',TicksInds)
%     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
    set(gca,'XTickLabel',dates)
    if k == 1
        title('\beta_t "real-time" estimates')
    end
    ymax = 3;
    ylim([0.5 ymax])
    hold on
    yis1 = ymax/100:ymax/100:ymax*2/5;
    yis2 = ymax*3/5:ymax/100:ymax;

%     try
%         plot(Data.Instants(8-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(8-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(8-1),ymax*1/2,'sch. closure','HorizontalAlignment','center')
%         plot(Data.Instants(14-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(14-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(14-1),ymax*1/2,'end of holidays','HorizontalAlignment','center')
%     end
    hold off

    
    
%     subplot(4,3,3*(k-1)+3)
%     
%     try
%         temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
%     catch
%         temp = diag(Res.Thetas(Parameters.gammam1.Index,:))*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
%     end
%     ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
%     hold on
%     ciplot(quantile(temp,0.25),quantile(temp,0.75),[100,153,251]/255)
%     plot(mean(temp),'k','LineWidth',2)
%     t = 0:21:567;
%     plot(t,1*ones(size(t)),'--k','LineWidth',2)
%     try
%         plot(Res.RtPath,'g','LineWidth',2) 
%     end
%     hold off
%     xlim([0 568])
%     TicksInds = [21*6 21*13 21*20 21*27]+1;
%     set(gca,'XTick',TicksInds)
% %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
% %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%     set(gca,'XTickLabel',dates)
%     title('R_t')
%     ymax = 3;
%     ylim([0.5 ymax])
%     hold on
%     yis1 = ymax/100:ymax/100:ymax*2/5;
%     yis2 = ymax*3/5:ymax/100:ymax;
%     try
%         plot(Data.Instants(8-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(8-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(8-1),ymax*1/2,'sch. closure','HorizontalAlignment','center')
%         plot(Data.Instants(14-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(14-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(14-1),ymax*1/2,'end of holidays','HorizontalAlignment','center')
%     end
%     hold off
end
    
    
    %%%%%%%%%%%ResMIF
%     Res = RessMIF{k};
%     Parameters = Res.ParamsMIF;
%     
%     Paths = Res.Paths;
%     
%     dates = {};
%     delta = floor(length(Res.Data.Dates)/Resol);
%     Resol = floor(length(Res.Data.Dates)/delta);
%     inds = delta:delta:delta*Resol;
%     for i = 1:length(inds)
%         dates{i} = [num2str(Res.Data.Dates{inds(i)}.Day) ' ' Res.Data.Dates{inds(i)}.MonthInLetters];
%     end
%    
%     
%     subplot(4,3,3*(k-1)+3)
%     
%     try
%         temp = diag(Res.Thetas(Parameters.gamma.Index,:).^-1)*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
%     catch
%         temp = diag(Res.Thetas(Parameters.gammam1.Index,:))*squeeze(exp(Paths(:,6,:)).*Paths(:,1,:))/Res.Parameters.TotalPopulation;
%     end
%     temp = squeeze(exp(Paths(:,6,:)));
%     ciplot(quantile(temp,0.025),quantile(temp,0.975),[172,215,255]/255)
%     hold on
%     ciplot(quantile(temp,0.25),quantile(temp,0.75),[100,153,251]/255)
%     plot(mean(temp),'k','LineWidth',2)
%     t = Data.Instants-Data.NbComputingSteps(end);
%     plot(t,1*ones(size(Data.Instants)),'--k','LineWidth',2)
%     try
%         plot(Res.RtPath,'g','LineWidth',2) 
%     end
%     hold off
%     xlim([0 567])
%     TicksInds = [delta*21*6 delta*21*13 delta*21*20];
%     set(gca,'XTick',TicksInds)
%     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
%     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%     set(gca,'XTickLabel',dates)
%     title('\beta_t')
%     ymax = 3;
%     ylim([0.5 ymax])
%     hold on
%     yis1 = ymax/100:ymax/100:ymax*2/5;
%     yis2 = ymax*3/5:ymax/100:ymax;
%     try
%         plot(Data.Instants(8-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(8-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(8-1),ymax*1/2,'sch. closure','HorizontalAlignment','center')
%         plot(Data.Instants(14-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(14-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(14-1),ymax*1/2,'end of holidays','HorizontalAlignment','center')
%     end
%     hold off
    
%     ciplot(quantile(squeeze(exp(Paths(:,6,:))),0.025),quantile(squeeze(exp(Paths(:,6,:))),0.975),[172,215,255]/255)
%     hold on
%     ciplot(quantile(squeeze(exp(Paths(:,6,:))),0.25),quantile(squeeze(exp(Paths(:,6,:))),0.75),[100,153,251]/255)
%     plot(mean(squeeze(exp(Paths(:,6,:)))),'k','LineWidth',2)
%     t = Data.Instants-Data.NbComputingSteps(end)*Res.Parameters.ComputationTStep;
%     try
%         plot(Res.RtPath,'g','LineWidth',2) 
%     end
%     hold off
%     xlim([0 1400])
%     TicksInds = [420 910 1400];
%     set(gca,'XTick',TicksInds)
% %     set(gca,'XTick',[t(end)/Resol:t(end)/Resol:t(end)])
% %     set(gca,'XTickLabel',['jun';'jul';'aug';'sep';'oct';'nov';'dec';'jan';'feb';'mar';'apr';'may'])
%     set(gca,'XTickLabel',dates)
%     title('\beta_t')
%     ymax = 3;
%     ylim([0.6 ymax])
%     hold on
%     yis1 = ymax/100:ymax/100:ymax*2/5;
%     yis2 = ymax*3/5:ymax/100:ymax;
%     try
%         plot(Data.Instants(8-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(8-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(8-1),ymax*1/2,'sch. closure','HorizontalAlignment','center')
%         plot(Data.Instants(14-1)*ones(size(yis1)),yis1,'r')
%         plot(Data.Instants(14-1)*ones(size(yis2)),yis2,'r')
%         text(Data.Instants(14-1),ymax*1/2,'end of holidays','HorizontalAlignment','center')
%     end
%     hold off
end


Parameters = struct();



SavePath = 'S:\Results\';
Res = load([SavePath 'Marc_ForPresEst_WithPaths.mat']);
Res = Res.Res3;
tmp = mean(Res.Thetas');
Parameters = Res.Parameters;
Parameters.Problem = 'MarcFlu';
Parameters.NbVariables = 7;
Parameters.SigmaObs = 0.1;
Parameters.DiffusionType = 'IBM';
Parameters.ObservationLength = 7*35;
Parameters.ComputationTStep = 0.1;
Parameters.TotalPopulation = 100000;
Parameters.betaderinit.Value = 0;
Parameters.betaderinit.Min = -10^14;
Parameters.betaderinit.Max = 10^14;
Parameters.betaderinit.MinLim = -100;
Parameters.betaderinit.MaxLim =  100;
Parameters.betaderinit.Estimated = 0;
Parameters.betaderinit.TransfType = 'Logit';
Parameters.betaderinit.Init = 1;
Parameters.betainit.Sample = 0;
Names = Res.Parameters.Names.Estimated;
for i = 1:length(tmp)
    Parameters.(Names{i}).Value = tmp(Res.Parameters.(Names{i}).Index);
end
Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);
TellParsValues(Parameters);


tmp = 1:4;
lengths = ceil(35/5*tmp);

for i = 2:length(lengths)
    Data.Observations = zeros(7,lengths(i));
    Data.Observations(5,:) = Weigthed(1:lengths(i))*10;
    Data.Instants = [1:lengths(i)]*7/Parameters.ComputationTStep;
    Data.ObservedVariables = 5*ones(1,length(Data.Instants));
    Data.NbComputingSteps = [0 diff(Data.Instants)];

    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsRealTime_IBM_' num2str(lengths(i)) 'weeks.mat'];
    FullSEIRinference(Data,Parameters,Name,1)
end

for i = 1:length(lengths)
    
    SavePath = 'S:\Results\';
    Name = [SavePath 'TestsRealTime_IBM_' num2str(lengths(i)) 'weeks.mat'];
    Res = load(Name);
    Res = Res.Res3;
    PlotMarc(Res,8)
    pause()
end













