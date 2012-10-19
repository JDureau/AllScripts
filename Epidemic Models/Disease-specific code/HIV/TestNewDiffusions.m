%% Testing New diffusions


% Gompertz
beta0 = 0.02;
delta = 0.5;
tinfl = 400;

x0 = log(beta0);
mu = beta0+delta;
k =  1/tinfl * log(log(mu)-x0);

tis = 0:0.01:1000;
xis = log(mu)+ (x0 - log(mu))*exp(-k*tis);
betas = exp(xis);

plot(tis,betas)


% Logistic
beta0 = 0.02;
delta = 0.5;
tinfl = 500;
steepness = 0.03;

mu = beta0 + delta;
x0 = 1/steepness*log(beta0/(mu-beta0));
k = -x0/(steepness*tinfl);

tis = 0:0.01:600;
xis = x0 + steepness*k*tis;
betas = mu * exp(xis)./(1+exp(xis));

clf
plot(tis,betas)
ylim([0 1])


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan';

cd('/Users/dureaujoseph/Dropbox/AllScriptsGit/')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])

SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan';

load([SavePath '/ParametersMysore.mat']);
Parameters.Problem = 'ImperialHIV';

Parameters.Problem = 'ImperialHIV';
Parameters.SigmaObs = 0.1;
Parameters.ComputationTStep = 0.5;
Parameters.ObservationLength = 25*12;
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
HIVModel.EKF_projection = @HIV_EKF_projection;
HIVModel.InitializeParameters = @HIV_Initialize;
HIVModel.SMC_projection = @HIV_SMC_projection;
HIVModel.LikFunction = 'normpdf(Variables(:,Data.ObservedVariables(:,IndTime)),Data.Observations(Data.ObservedVariables(:,IndTime),IndTime),sqrt(Data.Observations(Data.ObservedVariables(:,IndTime),IndTime)*(100-Data.Observations(Data.ObservedVariables(:,IndTime),IndTime))/400)).*(Res.WentOutOrNot)';


Mode = 1;
% 1 = Normal
% 2 = Biased
% 3 = Restricted
% 4 = AddingObs

ModelType = 'Bert';
% ModelType = 'Step';
% ModelType = 'Affine';
% ModelType = 'Sigm';
% ModelType = 'Logist';


NbTests = 400;

ampls = [];
inflpts = [];
slopes = [];amplints = [];
Ress = {};
cpt = 0;
ParsRecord = [];
ResGens = {};
NSamples = 100;
MaxValue = 0.6;
step = MaxValue/10;
Boxlimits = 0:step:MaxValue;
MaxEff = ceil(NSamples/(length(Boxlimits)-1));
BoxCpt = zeros(1,length(Boxlimits)-1);
while cpt < NSamples;
    try
        
        
        tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
        Parameters.NbTSteps = length(tis);

        
        if strcmp(ModelType,'Bert')
            asympt = 0.3+0.59*rand(1,1)^(((100-cpt)/100)^2);
%             m = 1+30*rand(1,1)^0.8; % when m not informed
            m = 6+60*rand(1,1);%^(((100-cpt)/100)^2);
            baseline = min(0.3*rand(1,1)^((2*(cpt)/100)^2),0.999*asympt*m^(1/(1-m)));
%             baseline = min(0.1*rand(1,1),0.999*asympt*m^(1/(1-m)));

            % m = 0.9*rand(1,1);
            inflpt = (300 + 280*rand(1,1)^(((100-cpt)/100)^2))*Parameters.ComputationTStep;

            B  = (1 - (baseline/asympt)^(1-m));
            k = (1/inflpt*log(B/(1-m)));

            Parameters.DiffusionType = 'Bertallanfy';


    %         xis =
    %         Parameters.ComputationTStep:Parameters.ComputationTStep:Param
    %         eters.ObservationLength;

            xis = [];
            xis(1) = (baseline^(1-m)-asympt^(1-m))/(1-m);
            sigma = 0;
            for i = 2:length(tis)
                xis(i) = xis(i-1) -k*xis(i-1)*Parameters.ComputationTStep + sigma*xis(i-1)*sqrt(Parameters.ComputationTStep).*randn(1,1);
            end
            Fts =  ((1-m)*xis+asympt^(1-m)).^(1/(1-m));
            plot(Fts)
            
        elseif strcmp(ModelType,'Logist')
            Parameters.DiffusionType = 'Bertallanfy';

            asympt = 0.4+0.59*rand(1,1);
            baseline = min(0.6*rand(1,1),0.999*asympt*m^(1/(1-m)));
            r = 0.05*rand(1,1);%*^1/10;
            
            Fts =  asympt*baseline*exp(r*tis)./(asympt+baseline*(exp(r*tis*Parameters.ComputationTStep)-1));
        elseif strcmp(ModelType,'Step')

            asympt = 0.4+0.59*rand(1,1);
%             baseline = min(0.6*rand(1,1),asympt);
            baseline = min(0.1*rand(1,1),asympt);
            inflpt = (300 + 280*rand(1,1)^(0.8))*Parameters.ComputationTStep;   
            
            Fts =  baseline + (asympt-baseline).*(tis>inflpt);
          
        elseif strcmp(ModelType,'Affine')

            asympt = 0.4+0.59*rand(1,1);
            baseline = min(0.6*rand(1,1),asympt);
            
            Fts =  baseline + (asympt-baseline).*(tis-tis(1))./(tis(end)-tis(1));
           
          
        elseif strcmp(ModelType,'Sigm')
        
            ampl = rand(1,1);
            baseline = rand(1,1)*(1-ampl);
            asympt = baseline+ampl;
            inflpt = (300 + 250*rand(1,1))*Parameters.ComputationTStep;
            steepness = rand(1,1)*50;

            Fts = baseline+ampl*Sigmoid((tis-inflpt)/steepness);
        end
        
        ampls(cpt+1) = Fts(end)-Fts(1);
        tmp = diff(Fts);
        [b,ind] = max(tmp);
        inflpts(cpt+1) = ind;
%         slopes(cpt+1) = asympt*k*m^(m/(1-m));
        amplints(cpt+1) = Fts(end)-Fts(456);
        
        plot(Fts)
        ylim([0 1])
        hold on
        plot([inflpts(cpt+1) inflpts(cpt+1)],[0 1])
        [b,ind]= max(diff(Fts));
%         plot([ind ind],[0 1],'k')
        hold off
%         title(r)
        
        
        if strcmp(ModelType,'Bert')
            Parameters.BRmu.Value = asympt;
            Parameters.BRmm1.Value = m-1;
            Parameters.BRbase.Value = baseline;
            Parameters.BRtinfl.Value = inflpt;
        elseif or(strcmp(ModelType,'Logist'),strcmp(ModelType,'Sigm'))
             Parameters.BRmm1.Value = 1;
             Parameters.BRbase.Value = 0.1;
             Parameters.BRtinfl.Value = 100;
             Parameters.BRmu.Value = 0.8;
        end
        Parameters.MultNoise = 0.05;
        Parameters = DefineEstimatedParametersIndexes(Parameters);
        Parameters = DefineTransfFunctions(Parameters);
        Parameters = DefinePriors(Parameters);
        Parameters = UpdateParsNoTransfToTransf(Parameters);
        ResGen = GenerateDataForGivenSig(Data,Parameters,HIVModel,Fts);
%         if(Fts(end)-Fts(1)<0.3)
%             if rand(1,1)<0.7
%                 die
%             end
%         end
               pause(0.01)

        if strcmp(ModelType,'Bert')
            x = Fts(end)-Fts(456);
            indbox = find(and(x<Boxlimits(2:end),x>=Boxlimits(1:end-1)),1);
            if or(isempty(indbox),BoxCpt(indbox)>=MaxEff)
                die
            else
                indbox(1)
                BoxCpt(indbox) = BoxCpt(indbox)+1;
                disp(indbox)
                disp(BoxCpt)
                cpt = cpt+1;
                disp([cpt sum(BoxCpt)])
            end
        end
                
      
        
       pause(0.01)
       disp(cpt)

       ResGens{end+1} = ResGen;
    end
end
subplot(3,1,1)
hist(ampls)
subplot(3,1,2)
hist(inflpts)
subplot(3,1,3)
hist(amplints)
	
if strcmp(ModelType,'Bert')
    save([SavePath '/ResGenBR_v2.mat'],'ResGens')
elseif strcmp(ModelType,'Logist')
    save([SavePath '/ResGenLogist.mat'],'ResGens')
elseif strcmp(ModelType,'Sigm')
    save([SavePath '/ResGenSigm.mat'],'ResGens')
elseif strcmp(ModelType,'Step')
    save([SavePath '/ResGenStep.mat'],'ResGens')
elseif strcmp(ModelType,'Affine')
    save([SavePath '/ResGenAffine.mat'],'ResGens')
end


tmp = [];
for i = 1:100
    tmp(i) = ResGens{i}.Data.Fts(1);
end

for i = 1:100
    plot(ResGens{i}.Data.Fts)
    ylim([0 1])
    pause()
end


load([SavePath '/ResGenBR.mat'])



Data = ResGens{87}.Data;
Parameters = ResGens{87}.Parameters;
Parameters.DiffusionType = 'Bertallanfy';


Parameters.BRm.Value = 3;
Parameters.BRm.Min = -10^14;
Parameters.BRm.Max = 10^14;
Parameters.BRm.MaxLim = 50;
Parameters.BRm.MinLim = 1.01;
Parameters.BRm.Estimated = 1;
Parameters.BRm.TransfType = 'Logit';
Parameters.BRm.Init = 1;

Parameters.BRmu.Value = 0.8;
Parameters.BRmu.Min = -10^14;
Parameters.BRmu.Max = 10^14;
Parameters.BRmu.MaxLim = 0.99;
Parameters.BRmu.MinLim = 0.4;
Parameters.BRmu.Estimated = 1;
Parameters.BRmu.TransfType = 'Logit';
Parameters.BRmu.Init = 1;

Parameters.BRbase.Value = 0.2;
Parameters.BRbase.Min = -10^14;
Parameters.BRbase.Max = 10^14;
Parameters.BRbase.MaxLim = 0.4;
Parameters.BRbase.MinLim = 0.01;
Parameters.BRbase.Estimated = 1;
Parameters.BRbase.TransfType = 'Logit';
Parameters.BRbase.Init = 1;

Parameters.BRtinfl.Value = 400;
Parameters.BRtinfl.Min = -10^14;
Parameters.BRtinfl.Max = 10^14;
Parameters.BRtinfl.MinLim = 1;
Parameters.BRtinfl.MaxLim = 600;
Parameters.BRtinfl.Estimated = 1;
Parameters.BRtinfl.TransfType = 'Logit';
Parameters.BRtinfl.Init = 1;

Parameters = DefineEstimatedParametersIndexes(Parameters);
Parameters = DefineTransfFunctions(Parameters);
Parameters = DefinePriors(Parameters);
Parameters = UpdateParsNoTransfToTransf(Parameters);

Names = Parameters.Names.Estimated;
NamesEst = {};
for i = 1:length(Names)
    if not(or(strcmp(Names{i},'SigmaRW'),strcmp(Names{i},'InitialFt')))
        NamesEst{end+1} = Names{i};
    end
end
Initialization = [];
for i = 1:length(NamesEst)
    Initialization(i) = Parameters.(NamesEst{i}).TransfValue ;
end
[x,fval,exitflag,output] = fminsearch(@(x) SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel),Initialization,optimset('MaxIter',10000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',20000));
for i = 1:length(NamesEst)
    Parameters.(NamesEst{i}).TransfValue = (x(i));
end
finalx = x;
Parameters = UpdateParsTransfToNoTransf(Parameters);
TellParsValues(Parameters)
InitParameters = Parameters;

fun =  @(x)SquarredErrorBRtoMinimize(x,NamesEst,Data,Parameters,HIVModel);
% [HD,err] = hessian(fun,Parameters.Finalx);
% Hess = -HD;
[HD,err] = hessdiag(fun,finalx);       
Hess = diag(-HD);

NbIts = 1000;
Parameters = InitParameters;
Cov = 2.38^2/length(finalx)*(-Hess)^-1;
Parameters.NamesEst = NamesEst;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts)


Parameters = Res.Parameters;
Cov = 2.38^2/length(finalx)*cov(Res.xs');
Parameters.NamesEst = NamesEst;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,NbIts)



Parameters = Res.Parameters;
Cov = 2.38^2/length(finalx)*cov(Res.xs');
Parameters.NamesEst = NamesEst;
Res = RunMCMCdeterm(Parameters,Data,HIVModel,Cov,2000)

PlotResHIV(Res,Res.Parameters)


SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan';

load([SavePath '/VIH_PlayingWithSigm_ResultsBertallanfy__1_1.mat'])


load([SavePath '/VIH_PlayingWithSigmStep10_ResultsAdd__1_1.mat'])
Res = Ress{1}
Res.Parameters.TypeWork = 'Boston Examples';
Res.Parameters.DiffusionType='Add';
PlotResHIV(Res,Res.Parameters)
% 
% figure(1)
load([SavePath '/PostThreeMethods.mat' ])




PlotPostForq('ROC',0.5,'Delta2003',0.2)





figure(1)
PlotPostForq('Mean')
% PlotPostForq(Res,'ROC',0.5)

load([SavePath '/PostThreeMethodsAffine.mat' ])
figure(2)
PlotPostForq(Res,'Mean',0.5)
% PlotPostForq(Res,'ROC',0.5)

load([SavePath '/PostThreeMethodsStep.mat' ])
figure(3)
PlotPostForq(Res,'Mean',0.5)
% PlotPostForq(Res,'ROC',0.5)

load([SavePath '/PostThreeMethods.mat' ])
figure(4)
PlotPostForq(Res,'Mean',0.5)
% PlotPostForq(Res,'ROC',0.5)

load([SavePath '/PostThreeMethodsAffine.mat' ])
figure(5)
PlotPostForq(Res,'Mean',0.5)
% PlotPostForq(Res,'ROC',0.5)

load([SavePath '/PostThreeMethodsStep.mat' ])
figure(6)
PlotPostForq(Res,'Mean',0.5)
% PlotPostForq(Res,'ROC',0.5)


load([SavePath '/PostDetmInf.mat' ])

figure(2)
load([SavePath '/PostThreeMethodsConstr.mat' ])

load([SavePath '/PostThreeMethodsLogist.mat' ])

figure(1)
Mode = 'ROC';

figure(2)
Mode = 'Mean';

ResAggr = Res;

figure(1)
PlotPostForq(Res,'Mean',0.5)
figure(2)
PlotPostForq(Res,'ROC',0.5)

PlotPostForq('ROC',0.5,'',0.2)

ToBar = [];
xis = [];
TempRes = Res;
tab = TempRes.partampls2003;
for  i = 1:6
    q1 = quantile(tab(4,:),0.16*(i-1));
    q2 = quantile(tab(4,:),0.16*(i));
    inds = find(and((tab(4,:)<q2),(tab(4,:)>q1)));
    ToBar(i,1) = mean(tab(1,inds) - tab(4,inds));
    ToBar(i,2) = mean(tab(3,inds) - tab(4,inds));
    ToBar(i,3) = mean(tab(2,inds) - tab(4,inds));
    xis(i) = mean([q1 q2]);
end
bar(1:5,ToBar(2:end,:)./repmat(xis(2:end)',1,3))
% bar(xis,ToBar)
xlabel('"True" \DeltaCU')
ylabel('Relative bias')
legend('dBR','sBR','BM')
set(gca,'XTick',[0:(600)/5:(600)])
set(gca,'XTickLabel',['1985';'1990';'1995';'2000';'2005';'2010'])

Title('Relative bias in CU estimate as a function of \DeltaCU')

figure(7)

clf
for i = 1:length(ResAggr.ampls)
%     plot(ResGens{i}.Data.Fts)
%     hold on
%     plot(2*[ResGens{i}.Parameters.BRtinfl.Value ResGens{i}.Parameters.BRtinfl.Value],[0 1],'g')
%     hold off
    

    plot(squeeze(ResAggr.sigms(4,i,:)),'g')
    hold on
    plot(squeeze(ResAggr.sigms(1,i,:)),'k')
    plot(squeeze(ResAggr.sigms(3,i,:)),'b')
    plot(squeeze(ResAggr.sigms(2,i,:)),'r')
    plot(squeeze(ResAggr.sigms(5,i,:)),'y')
    plot(squeeze(ResAggr.sigms(6,i,:)),'c')
    plot([ResAggr.ParsMeds(1,16,i) ResAggr.ParsMeds(1,16,i)],[0 1],'k')
    plot([ResAggr.ParsMeds(3,16,i) ResAggr.ParsMeds(3,16,i)],[0 1],'b')
%     plot(2*[ResGens{i}.Parameters.BRtinfl.Value ResGens{i}.Parameters.BRtinfl.Value],[0 1],'g')
    hold off
%     title([ResAggr.ParsMeds(1,16,i) ResAggr.ParsMeds(3,16,i) ResAggr.ParsMeds(4,16,i)] )
    ylim([0 1])
%     legend('true','det','Ber','Add')
    pause()
end



PlotPostForAll(Res,'Mean')





q = 0.5;
Stats = {'Ampl. after 1995','Ampl after 2003','2010 value'};
Tabs = {};



% tmpinds = find(not(isnan(Res.ampls(1,:))));
% Res.ampls = Res.ampls(:,tmpinds);
% Res.partampls1995 = Res.partampls1995(:,tmpinds);
% Res.partampls2003 = Res.partampls2003(:,tmpinds);
% Res.asympts = Res.asympts(:,tmpinds);
% tmp = find(not(Res.NbSamples(:,1)==0));
% Res.NbSamples = Res.NbSamples(tmp,:);
% inds = find(and(Res.NbSamples(:,2)>50, Res.NbSamples(:,3)>20));
% length(inds)
inds = 1:size(Res.ampls,2);


if strcmp(Mode,'ROC')
    Tabs{1} = Res.amplsROC(:,inds);
    Tabs{2} = Res.partampls1995ROC(:,inds);
    Tabs{3} = Res.partampls2003ROC(:,inds);
    Tabs{4} = Res.asymptsROC(:,inds);
elseif strcmp(Mode,'Post')
    Tabs{1} = Res.amplsPost(:,inds);
    Tabs{2} = Res.partampls1995Post(:,inds);
    Tabs{3} = Res.partampls2003Post(:,inds);
    Tabs{4} = Res.asymptsPost(:,inds);
else
    Tabs{1} = Res.ampls(:,inds);
    Tabs{2} = Res.partampls1995(:,inds);
    Tabs{3} = Res.partampls2003(:,inds);
    Tabs{4} = Res.asympts(:,inds);
end





Biases = [];
MSEs = [];
Vars = [];
ROCs = [];
inds = [1 2];
Pts = {};
for i = 1:4
    for j = [1 2]%1:3
        Corrs(inds(j),i) = corr(Tabs{i}(4,:)',Tabs{i}(j,:)');
        Biases(inds(j),i) = mean(Tabs{i}(j,:)-Tabs{i}(4,:)); 
        MSEs(inds(j),i) = mean((Tabs{i}(4,:)-Tabs{i}(j,:)).^2);
        Vars(inds(j),i) = MSEs(inds(j),i)-Biases(inds(j),i).^2;
        
        quan = quantile(Tabs{i}(4,:),q);
%         med = median(Tabs{i}(4,:));
        indsinf = find(Tabs{i}(4,:)<=quan);
        indssup = find(Tabs{i}(4,:)>quan);
        xs = [];
        ys = [];
        all = [Tabs{i}(1,:) Tabs{i}(2,:) Tabs{i}(3,:)];
        m = min(all);
        M = max(all);
%         deltap = sort((all(1:end-1)+all(2:end))/2);
        deltap = sort(all);
        deltap = (deltap(1:end-1)+deltap(2:end))/2;
        deltap = [(min(all)-mean(diff(deltap))) deltap (max(all)+mean(diff(deltap)))];
        for k = 1:length(deltap)
            xs(k) = sum(Tabs{i}(j,indsinf)>deltap(k))/length(indsinf);
            ys(k) = sum(Tabs{i}(j,indssup)>deltap(k))/length(indssup);
        end
        ROCs(inds(j),i) = 0;
        for k = 1:length(deltap)-1
            ROCs(inds(j),i) = ROCs(inds(j),i) + (xs(k)-xs(k+1))*(ys(k+1)+ys(k))/2;
        end
        ROCs(inds(j),i) = ROCs(inds(j),i)-0.5;
%         ROCs(inds(j),i) = -diff(xs)*(ys(1:end-1)'+ys(2:end)')/2-0.5;
        Pts{i,inds(j)}=[xs;ys];
    end
end

s = 16

subplot(4,3,1:3)
bar(Biases(:,2:4)')
set(gca,'XTickLabel',Stats,'FontWeight','bold','FontSize',s)
ylim([-0.2 0.2])
ylabel('Bias')

subplot(4,3,4:6)
bar(Vars(:,2:4)')
set(gca,'XTickLabel',Stats,'FontWeight','bold','FontSize',s)
ylabel('Variance')

subplot(4,3,7:9)
bar(MSEs(:,2:4)')
ylim([0 0.04])
set(gca,'XTickLabel',Stats,'FontWeight','bold','FontSize',s)
ylim([0 0.04])
ylabel('MSE')
legend('Parametric','Diffusion')

for i = 2:4
   subplot(4,3,8+i) 
   plot(Pts{i,1}(1,:),Pts{i,1}(2,:),'b')
   hold on
   plot(Pts{i,2}(1,:),Pts{i,2}(2,:),'r')
   hold off
   xlabel('False pos. rate','FontWeight','bold','FontSize',s)
   ylabel('True pos. rate','FontWeight','bold','FontSize',s)
end





axes('Position',[0 0 1 1],'Visible','off');
text(0.45,0.195,'Title')

clf
for i = 1:length(ResAggr.ampls)
    plot(squeeze(ResAggr.sigms(4,i,:)),'g')
    hold on
    plot(squeeze(ResAggr.sigms(1,i,:)),'k')
    plot(squeeze(ResAggr.sigms(3,i,:)),'b')
    plot(squeeze(ResAggr.sigms(2,i,:)),'r')
    hold off
    title(i)
    ylim([0 1])
%     legend('true','det','Ber','Add')
    pause()
end

indsinfl = [];
for i = 1:length(ResAggr.ampls)
    [b,indsinfl(i)] = max(diff(squeeze(ResAggr.sigms(1,i,:))));
end
hist(indsinfl)
subplot(3,1,1)
plot(ResAggr.ampls(4,:),ResAggr.ampls(1,:),'.')
title(corr(ResAggr.ampls(4,:)',ResAggr.ampls(1,:)'))
subplot(3,1,2)
plot(ResAggr.ampls(4,:),ResAggr.ampls(3,:),'.')
title(corr(ResAggr.ampls(4,:)',ResAggr.ampls(3,:)'))
subplot(3,1,3)
plot(ResAggr.ampls(4,:),ResAggr.ampls(2,:),'.')
title(corr(ResAggr.ampls(4,:)',ResAggr.ampls(2,:)'))



for i = 1:5
    try
        load([SavePath '/ResGenBRBert_BRdeterm_' num2str(i) '.mat'])
        ResDet = Res;
        load([SavePath '/VIH_PlayingWithSigm_ResultsBertallanfy__1_' num2str(i) '.mat'])
        ResBer = Ress{1};
        ResBer.Parameters.TypeWork = 'Boston Examples';
        load([SavePath '/VIH_PlayingWithSigm_ResultsAdd__1_' num2str(i) '.mat'])
        ResAdd = Ress{1};
        ResAdd.Parameters.TypeWork = 'Boston Examples';
        load([SavePath '/ResGenBR_BRdetermConstr_' num2str(i) '.mat'])
        ResDetConstr = Res;
        load([SavePath '/VIH_PlayingWithSigm_ResultsBertallanfyConstr__1_' num2str(i) '.mat'])
        ResBerConstr = Ress{1};
        ResBerConstr.Parameters.TypeWork = 'Boston Examples';
        load([SavePath '/VIH_PlayingWithSigm_ResultsAddConstr__1_' num2str(i) '.mat'])
        ResAddConstr = Ress{1};
        ResAddConstr.Parameters.TypeWork = 'Boston Examples';
        
        if or(not(ResDet.Data.BuiltTraj(1:585,9) == ResAdd.Data.BuiltTraj(1:585,9)),not(ResDet.Data.BuiltTraj(1:585,9) == ResBer.Data.BuiltTraj(1:585,9)))
          disp('NOT OK')
          die
        end

%             Names = {};
%             for j = 2:13
%                 Names{end+1} = ResAdd.Parameters.Names.Estimated{j};
%             end
%             TempRes = ResBer;
%             Names = TempRes.Parameters.Names.Estimated;
%             
%             for j = 1:length(Names)
%                 subplot(3,6,j)
%                 plot(TempRes.Thetas(TempRes.Parameters.(Names{j}).Index,:))
% %                 plot(ResDet.xs(ResAdd.Parameters.(Names{j}).Index,:))
% 
%     %             plot(ResBer.Thetas(ResBer.Parameters.(Names{j}).Index,:),ResAdd.Paths(:,3,300),'.')
%     %             [fi,xi] = ksdensity(ResAdd.Thetas(ResAdd.Parameters.(Names{j}).Index,:));
%     %             plot(xi,fi,'r')
%     %             hold on
%     % % %             [fi,xi] = ksdensity(ResDet.Thetas(ResDet.Parameters.(Names{i}).Index,:));
%     % % %             plot(xi,fi,'b')
%     %             [fi,xi] = ksdensity(ResBer.Thetas(ResBer.Parameters.(Names{j}).Index,:));
%     %             plot(xi,fi,'g')
%     %             hold off
%                 title(Names{j})
% %                 v = ResDet.xs(ResAdd.Parameters.(Names{j}).Index,:);
%                 v = TempRes.Thetas(TempRes.Parameters.(Names{j}).Index,:);
%                 tmp = Autocorrelation(v,100);
% %                 title(length(v)/(1+2*sum(tmp(2:end))))
%             end
        
%         subplot(2,1,1)
%         plot(ResBer.LogPosts(20:end))
%         subplot(2,1,2)
%         plot(ResBer.LogLiks(20:end))
%         subplot(3,1,3)
%         plot(ResDet.    )
        
        ResDetConstr.Parameters.PlotIndex = 1;
        PlotResHIV(ResDetConstr,ResDetConstr.Parameters)
        title('Deterministic Berthallanfy Constr','FontWeight','bold')

        ResBerConstr.Parameters.PlotIndex = 2;
        PlotResHIV(ResBerConstr,ResBerConstr.Parameters)
        title('Stochastic Berthallanfy Constr','FontWeight','bold')

        ResAddConstr.Parameters.PlotIndex = 3;
        PlotResHIV(ResAddConstr,ResAddConstr.Parameters)
        title('Brownian Motion Constr','FontWeight','bold')
        
        ResDet.Parameters.PlotIndex = 4;
        PlotResHIV(ResDet,ResDet.Parameters)
        title('Deterministic Berthallanfy','FontWeight','bold')

        ResBer.Parameters.PlotIndex = 5;
        PlotResHIV(ResBer,ResBer.Parameters)
        title('Stochastic Berthallanfy','FontWeight','bold')

        ResAdd.Parameters.PlotIndex = 6;
        PlotResHIV(ResAdd,ResAdd.Parameters)
        title('Brownian Motion','FontWeight','bold')
        
        

%         figure(4)
%         plot(squeeze(ResAggr.sigms(4,i,:)),'g')
%         hold on
%         plot(squeeze(ResAggr.sigms(1,i,:)),'k')
%         plot(squeeze(ResAggr.sigms(3,i,:)),'b')
%         plot(squeeze(ResAggr.sigms(2,i,:)),'r')
%         hold off
%         ylim([0 1])
%         
        ToSave = '/Users/dureaujoseph/Dropbox/Taf/Notes/Avahan';
        saveas(gcf,[ToSave '/ExampleThreeMethods' num2str(i)], 'jpg')
        i
        pause()
    end
end



%% Real Data




cd('/Users/dureaujoseph/AllScripts')
addpath([pwd '/General Tools'])
addpath([pwd '/Toolboxes'])
addpath([pwd '/Epidemic Models/Generic PMCMC tools'])
addpath([pwd '/Epidemic Models/Disease-specific code'])
addpath([pwd '/Epidemic Models/Disease-specific code/HIV'])
SavePath = '/Users/dureaujoseph/Documents/PhD_Data/Avahan/';


DistrictNames = {'Mysore_3rounds','Belgaum_3rounds','Bellary','EastGodavry','Guntur','Hyderabad','Shimoga'};

for i = 1:length(DistrictNames)
    clf
    load([SavePath '/HIV_' DistrictNames{i} '_d2.mat'])
    Res.Parameters.PlotIndex = 1;
    Res.Parameters.TypeWork='Normal';
    PlotResHIV(Res,Res.Parameters)
%     title([DistrictNames{i} ' BM'],'FontWeight','bold')
%     load([SavePath '/HIV_' DistrictNames{i} 'Bertallanfy.mat'])
%     Res.Parameters.PlotIndex = 2;
%     Res.Parameters.TypeWork='Boston Examples';
%     PlotResHIV(Res,Res.Parameters)
%     title([DistrictNames{i} ' sBR'],'FontWeight','bold')
    DistrictNames{i}
    pause()
end



load([SavePath '/HIV_Mysore_3roundsDetBertallanfy.mat'])
ResDet = Res;
load([SavePath '/HIV_Mysore_3roundsAdd.mat'])
ResAdd = Res;
load([SavePath '/HIV_Mysore_3roundsBertallanfy.mat'])
ResBer = Res;

Res.Parameters.TypeWork = 'Normal';
PlotResHIV(Res,Res.Parameters)

ResDet.Parameters.ObsMin =  ResBer.Parameters.ObsMin ;
ResDet.Parameters.ObsMax =  ResBer.Parameters.ObsMax ;


FtDet = mean(squeeze(ResDet.Paths(:,3,1:582)));
FtsDet = (squeeze(ResDet.Paths(:,3,1:582)));

tmp = squeeze(ResAdd.Paths(:,3,1:582));
FtAdd = mean(exp(tmp)./(1+exp(tmp)));
FtsAdd = (exp(tmp)./(1+exp(tmp)));

tmp = squeeze(ResBer.Paths(:,3,1:582));
beta0s = squeeze(ResBer.Thetas(ResBer.Parameters.BRbase.Index,:));
mus = squeeze(ResBer.Thetas(ResBer.Parameters.BRmu.Index,:));
ms = squeeze(ResBer.Thetas(ResBer.Parameters.BRmm1.Index,:))+1;
ts = squeeze(ResBer.Thetas(ResBer.Parameters.BRtinfl.Index,:));
Bs = 1-(beta0s./mus).^(1-ms);
ks = 1./ts.*log(Bs./(1-ms));
FtBer = mean((repmat(1-ms',1,582).*tmp+repmat(mus.^(1-ms)',1,582)).^(repmat(1./(1-ms)',1,582)));
FtsBer = ((repmat(1-ms',1,582).*tmp+repmat(mus.^(1-ms)',1,582)).^(repmat(1./(1-ms)',1,582)));


[fi,xi] = ksdensity(FtsBer(:,582)-FtsBer(:,585-168));
plot(xi,fi,'g')
hold on
[fi,xi] = ksdensity(FtsAdd(:,582)-FtsAdd(:,585-168));
plot(xi,fi,'r')
[fi,xi] = ksdensity(FtsDet(:,582)-FtsDet(:,585-168));
plot(xi,fi,'b')
hold off

Names = {'Mysore_3rounds','Belgaum_3rounds','Bellary','Yevatmal','Guntur','Hyderabad','EastGodavry','Shimoga'};
ests = [];
randinds = 1:2000;
for i = 1:length(Names)
    disp(Names{i})
%     load([SavePath '/HIV_' Names{i} '_dBR.mat'])
    ResDet = Res;
    randinds = randsample(size(ResDet.Paths,1),min(size(ResDet.Paths,1),4000));
    ResDet.Paths = ResDet.Paths(randinds,:);
    ResDet.Thetas = ResDet.Thetas(:,randinds);
    load([SavePath '/HIV_' Names{i} '_Sigm.mat'])
    ResSigm = Res;
    randinds = randsample(size(ResSigm.Paths,1),min(size(ResSigm.Paths,1),4000));
    ResSigm.Paths = ResSigm.Paths(randinds,:);
    ResSigm.Thetas = ResSigm.Thetas(:,randinds);
    load([SavePath '/HIV_' Names{i} '.mat'])
    ResAdd = Res;
    randinds = randsample(size(ResAdd.Paths,1),min(size(ResAdd.Paths,1),4000));
    ResAdd.Paths = ResAdd.Paths(randinds,:,:);
    ResAdd.Thetas = ResAdd.Thetas(:,randinds);
%     load([SavePath '/HIV_' Names{i} 'Bertallanfy.mat'])
    ResBer = Res;
    randinds = randsample(size(ResBer.Paths,1),min(size(ResBer.Paths,1),4000));
    ResBer.Paths = ResBer.Paths(randinds,:,:);
    ResBer.Thetas = ResBer.Thetas(:,randinds);
%     load([SavePath '/HIV_' Names{i} 'Sigmoid.mat'])
    ResSigmSto = Res;
    randinds = randsample(size(ResSigmSto.Paths,1),min(size(ResSigmSto.Paths,1),4000));
    ResSigmSto.Paths = ResSigmSto.Paths(randinds,:,:);
    ResSigmSto.Thetas = ResSigmSto.Thetas(:,randinds);
    
%     
    indend = ResSigm.Data.Instants(end);
%     
%     FtDet = mean(squeeze(ResDet.Paths(:,1:indend)));
%     FtsDet = (squeeze(ResDet.Paths(:,1:indend)));
% 
    FtSigm = mean(squeeze(ResSigm.Paths(:,1:indend)));
    FtsSigm = (squeeze(ResSigm.Paths(:,1:indend)));
%     
    tmp = squeeze(ResAdd.Paths(:,3,1:indend));
    FtAdd = mean(exp(tmp)./(1+exp(tmp)));
    FtsAdd = (exp(tmp)./(1+exp(tmp)));
% 
%     tmp = squeeze(ResBer.Paths(:,3,1:indend));
%     beta0s = squeeze(ResBer.Thetas(ResBer.Parameters.BRbase.Index,:));
%     mus = squeeze(ResBer.Thetas(ResBer.Parameters.BRmu.Index,:));
%     ms = squeeze(ResBer.Thetas(ResBer.Parameters.BRmm1.Index,:))+1;
%     ts = squeeze(ResBer.Thetas(ResBer.Parameters.BRtinfl.Index,:));
%     Bs = 1-(beta0s./mus).^(1-ms);
%     ks = 1./ts.*log(Bs./(1-ms));
%     FtBer = mean((repmat(1-ms',1,indend).*tmp+repmat(mus.^(1-ms)',1,indend)).^(repmat(1./(1-ms)',1,indend)));
%     FtsBer = ((repmat(1-ms',1,indend).*tmp+repmat(mus.^(1-ms)',1,indend)).^(repmat(1./(1-ms)',1,indend)));
% 
%     tmp = squeeze(ResSigmSto.Paths(:,3,1:indend));
%     rate = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmrate.Index,:));
%     base = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmbase.Index,:));
%     mu = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmmu.Index,:));
%     tinfl = squeeze(ResSigmSto.Thetas(ResSigmSto.Parameters.Sigmtinfl.Index,:));;
%     c = 1./(1+exp(tinfl./rate));
%     b = (mu-base).*c;
%     a = base - b;
%     a = repmat(a',1,indend);
%     b = repmat(b',1,indend);
%     c = repmat(c',1,indend);
%     FtSigmSto = mean(a + b./(c.*(1+tmp)));
%     FtsSigmSto = (a + b./(c.*(1+tmp)));
    
%     q = 0.5;
%     disp('q50 delta')
%     disp(quantile(FtsDet(:,indend)-FtsDet(:,585-168),q))
%     disp(quantile(FtsBer(:,indend)-FtsBer(:,585-168),q))
%     disp(quantile(FtsAdd(:,indend)-FtsAdd(:,585-168),q))
%     disp(quantile(FtsSigm(:,indend)-FtsSigm(:,585-168),q))
%     disp(quantile(FtsSigmSto(:,indend)-FtsAdd(:,585-168),q))
%     ests(i,1,1) = quantile(FtsDet(:,indend)-FtsDet(:,585-168),q);
%     ests(i,1,2) = quantile(FtsBer(:,indend)-FtsBer(:,585-168),q);
%     ests(i,1,3) = quantile(FtsAdd(:,indend)-FtsAdd(:,585-168),q);
%     ests(i,1,4) = quantile(FtsSigm(:,indend)-FtsSigm(:,585-168),q);
%     ests(i,1,5) = quantile(FtsAdd(:,indend)-FtsAdd(:,585-168),q);
%     
%     q = 0.975;
%     disp('q5 delta')
%     disp(quantile(FtsDet(:,indend)-FtsDet(:,585-168),q))
%     disp(quantile(FtsBer(:,indend)-FtsBer(:,585-168),q))
%     disp(quantile(FtsAdd(:,indend)-FtsAdd(:,585-168),q))
%     disp(quantile(FtsSigm(:,indend)-FtsSigm(:,585-168),q))
%     disp(quantile(FtsSigmSto(:,indend)-FtsSigmSto(:,585-168),q))
%     ests(i,2,1) = quantile(FtsDet(:,indend)-FtsDet(:,585-168),q);
%     ests(i,2,2) = quantile(FtsBer(:,indend)-FtsBer(:,585-168),q);
%     ests(i,2,3) = quantile(FtsAdd(:,indend)-FtsAdd(:,585-168),q);
% 
%     disp(mean(FtsDet(:,indend)-FtsDet(:,585-168)))
%     disp(mean(FtsBer(:,indend)-FtsBer(:,585-168)))
%     disp(mean(FtsAdd(:,indend)-FtsAdd(:,585-168)))
%     disp(mean(FtsSigm(:,indend)-FtsSigm(:,585-168)))
%     disp(mean(FtsSigmSto(:,indend)-FtsSigmSto(:,585-168)))
%     
%     q = 0.5
%     disp('q50 2010')
%     quantile(FtsDet(:,indend),q)
%     quantile(FtsBer(:,indend),q)
%     quantile(FtsAdd(:,indend),q)
%     ests(i,3,1) = quantile(FtsDet(:,indend),q);
%     ests(i,3,2) = quantile(FtsBer(:,indend),q);
%     ests(i,3,3) = quantile(FtsAdd(:,indend),q);
%     
%     q = 0.05
%     disp('q5 2010')
%     quantile(FtsDet(:,indend),q)
%     quantile(FtsBer(:,indend),q)
%     quantile(FtsAdd(:,indend),q)
%     ests(i,4,1) = quantile(FtsDet(:,indend),q);
%     ests(i,4,2) = quantile(FtsBer(:,indend),q);
%     ests(i,4,3) = quantile(FtsAdd(:,indend),q);
%     
%     disp('mean delta')
%     ests(i,5,1) = mean(FtsDet(:,indend)-FtsDet(:,585-168));
%     ests(i,5,2) = mean(FtsBer(:,indend)-FtsBer(:,585-168));
%     ests(i,5,3) = mean(FtsAdd(:,indend)-FtsAdd(:,585-168));
%     
%     disp('mean 2010')
%     ests(i,6,1) = mean(FtsDet(:,indend));
%     ests(i,6,2) = mean(FtsBer(:,indend));
%     ests(i,6,3) = mean(FtsAdd(:,indend));
    
    
    clf
    ResAdd.Parameters.TypeWork = 'Boston Examples';
    ResAdd.Parameters.PlotIndex = 3;
    ResAdd.title = 'c) Estimated condom use (Brownian motion model)';

    PlotResHIV(ResAdd,ResAdd.Parameters)
%     ResDet.Parameters.TypeWork = 'Boston Examples';
%     ResDet.Parameters.PlotIndex = 2;
%     ResDet.title = 'dBR';
%     PlotResHIV(ResDet,ResDet.Parameters)
    ResSigm.Parameters.TypeWork = 'Boston Examples';
    ResSigm.Parameters.PlotIndex = 4;
    ResSigm.title = 'd) Estimated condom use (det. sigmoid model)';
    PlotResHIV(ResSigm,ResSigm.Parameters)
%     ResBer.Parameters.TypeWork = 'Boston Examples';
%     ResBer.Parameters.PlotIndex = 4;
%     PlotResHIV(ResBer,ResBer.Parameters)
%     ResSigmSto.Parameters.TypeWork = 'Boston Examples';
%     ResSigmSto.Parameters.PlotIndex = 2;
%     ResSigmSto.title=  'Stochastic sigmoid model';
%     PlotResHIV(ResSigmSto,ResSigmSto.Parameters)

    q = 0.5
    quantile(FtsSigm(:,min(indend,582))-FtsSigm(:,585-168),q)
    quantile(FtsAdd(:,min(indend,582))-FtsAdd(:,585-168),q)


%     pause()
end

save([SavePath '/AllRegionsEstimates.mat'],'ests')






% mean(FtsDet(:,582)-FtsDet(:,585-168))
% mean(FtsBer(:,582)-FtsBer(:,585-168))
mean(FtsSigm(:,582)-FtsSigm(:,585-168))
% mean(FtsSigmSto(:,582)-FtsSigmSto(:,585-168))
mean(FtsAdd(:,582)-FtsAdd(:,585-168))

q = 0.5
% quantile(FtsDet(:,582)-FtsDet(:,585-168),q)
% quantile(FtsBer(:,582)-FtsBer(:,585-168),q)
quantile(FtsSigm(:,582)-FtsSigm(:,585-168),q)
% quantile(FtsSigmSto(:,582)-FtsSigmSto(:,585-168),q)
quantile(FtsAdd(:,582)-FtsAdd(:,585-168),q)

mean(FtsDet(:,582))
mean(FtsBer(:,582))
mean(FtsAdd(:,582))

[fi,xi] = ksdensity(FtsDet(:,indend)-FtsDet(:,585-168));
plot(xi,fi,'g','LineWidth',2);
hold on
[fi,xi] = ksdensity(FtsBer(:,indend)-FtsBer(:,585-168));
plot(xi,fi,'r','LineWidth',2);
[fi,xi] = ksdensity(FtsAdd(:,indend)-FtsAdd(:,585-168));
plot(xi,fi,'LineWidth',2);
plot([median(FtsAdd(:,indend)-FtsAdd(:,585-168)) median(FtsAdd(:,indend)-FtsAdd(:,585-168))],[0 1.4],':','LineWidth',4)
plot([median(FtsDet(:,indend)-FtsDet(:,585-168)) median(FtsDet(:,indend)-FtsDet(:,585-168))],[0 1.4],':g','LineWidth',4)
plot([median(FtsBer(:,indend)-FtsBer(:,585-168)) median(FtsBer(:,indend)-FtsBer(:,585-168))],[0 1.4],':r','LineWidth',4)
hold off
legend('Det. BR','Sto. BR','Brownian motion')
ylabel('Posterior density','FontSize',14)
xlabel('Shift in CU during the intervention in Mysore','FontSize',14)
xlim([-1 1])
