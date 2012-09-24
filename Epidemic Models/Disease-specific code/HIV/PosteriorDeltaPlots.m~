
% Generate the plot for deltaCU priors


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



ampls = [];




ModelType = 'Bert';
sigs = [0 1];
pars = [5 5];

ModelType = 'Sigm';
sigs = [0 1];
pars = [5 5];

ModelType = 'BM';
sigs = [1];
pars = [10];

deltas = [];
NSamples = 5000;

ms = [];
for j = 1:length(sigs)
    cpt = 1;
    while cpt < NSamples
        try
            cpt

            tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
            Parameters.NbTSteps = length(tis);


            if strcmp(ModelType,'Bert')
                asympt = rand(1,1);
                m = rand(1,1)*200;%pars(j);%^(((100-cpt)/100)^2);

                baseline = rand(1,1)*0.999*asympt*m^(1/(1-m));
                if baseline > 0.999*asympt*m^(1/(1-m))
                    die
                end
                inflpt = 600*rand(1,1)*Parameters.ComputationTStep;

                B  = (1 - (baseline/asympt)^(1-m));
                k = (1/inflpt*log(B/(1-m)));
                ms(cpt) = m;

                xis = [];
                xis(1) = (baseline^(1-m)-asympt^(1-m))/(1-m);
                sigma = sigs(j)*rand(1,1)*0.5;
                for i = 2:length(tis)
                    xis(i) = xis(i-1) -k*xis(i-1)*Parameters.ComputationTStep + sigma*xis(i-1)*sqrt(Parameters.ComputationTStep).*randn(1,1);
                end
                Fts =  ((1-m)*xis+asympt^(1-m)).^(1/(1-m));


            elseif strcmp(ModelType,'Sigm')  
                rate = rand(1,1)*200;%pars(j);
                base = rand(1,1);
                mu = rand(1,1);
                tinfl = 600*rand(1,1)*Parameters.ComputationTStep;

                c = 1/(1+exp(tinfl/rate));
                b = (mu-base)*c;
                a = base - b;
                xis = [];
                xis(1) = exp(tinfl/rate);
                sigma = sigs(j)*rand(1,1)*0.5;
                for i = 2:length(tis)
                    xis(i) = xis(i-1) -1/rate*xis(i-1)*Parameters.ComputationTStep + sigma*xis(i-1)*sqrt(Parameters.ComputationTStep).*randn(1,1);
                end
                Fts = (a + b./(c.*(1+xis)));
            elseif strcmp(ModelType,'BM') 
                base = rand(1,1);
                sigma = rand(1,1)*0.5;%sigs(j);
                xis = [];
                xis(1) = base;
                for i = 2:length(tis)
                    xis(i) = xis(i-1) + sigma*sqrt(Parameters.ComputationTStep).*randn(1,1);
                end
                Fts = exp(xis)./(1+exp(xis));
            end
            deltas(j,cpt) = Fts(584)-Fts(416);
            if or(not(isreal(Fts(584)-Fts(416))),isnan(Fts(584)-Fts(416)))
                die
            end
            cpt = cpt+1;

        end
    end
end
save([SavePath '/PosteriorDeltas_' ModelType '.mat'],'deltas')



ModelType = 'Bert';
ModelType = 'Sigm';
ModelType = 'BM';
load([SavePath '/PosteriorDeltas_' ModelType '.mat'])

clf
dash = [0 1];
cols = [1  1];
for i = 1:size(deltas,1)
    [fi,xi] = ksdensity(deltas(i,:));
    if dash(i)
        if cols(i) == 1
            plot(xi,(fi),':k','LineWidth',1.5)
        else
            plot(xi,(fi),':','LineWidth',1.5)
        end
    else
        if cols(i) == 1
            plot(xi,(fi),'k','LineWidth',1.5)
        else
            plot(xi,(fi),'LineWidth',1.5)
        end
    end        
    
    inds = find(not(isnan(deltas(i,:))));
    
    q1 = quantile(deltas(i,inds),0.025);
    q2 = quantile(deltas(i,inds),0.5);
    q3 = quantile(deltas(i,inds),0.975);
    s = std(deltas(i,inds));
    disp([q1 q2 q3 s])
    hold on
    xlim([-0.5 0.5])
end
legend('\sigma = 0','\sigma < 0.5','m = 10, \sigma = 0','m = 10, \sigma = 0.2')
% legend('r = 10, \sigma = 0','r = 10, \sigma = 0.2','r = 100, \sigma = 0','r = 100, \sigma = 0.2')
% legend(' \sigma = 0.2','\sigma = 0.4')
hold off

subplot(3,1,1)
hist(ampls)
subplot(3,1,2)
hist(inflpts)
subplot(3,1,3)
hist(amplints)