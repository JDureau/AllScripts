function Pars = FitSigmoids(data,Parameters)


NbSamples = 1000;

ampls = [];
baselines = [];
endlines = [];
inflpts = [];
steepnesses = [];
sigmampls = [];
sigms  = [];

for i = 1:NbSamples
    ind = ceil(rand(1,1)*size(data,1));
    ind
    xinit(1) = data(ind,1); %baseline
    xinit(2) = max(data(ind,:))-min(data(ind,:)); % ampl
    [b,indinfl] = max(diff(data(ind,:)));
    xinit(3) = indinfl; % inflpt
    xinit(4) = b/Parameters.ComputationTStep;

    [x,fval,exitflag,output] = fminsearch(@(x) DistToSigm(x,data(ind,:),Parameters),xinit,optimset('MaxIter',50,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));

    xis = (1:size(data,2))*Parameters.ComputationTStep;
    sigm = x(1)+x(2)*Sigmoid((xis-x(3)*Parameters.ComputationTStep)/x(4));
        
    [b,indmed] = min(abs(sigm-(sigm(585)+sigm(1))/2));
    if indmed>200   
        
        ampls(end+1) = x(end)-x(1);
        baselines(end+1) = sigm(1);
        endlines(end+1) = sigm(min(length(sigm),583));
        sigmampls(end+1) = sigm(end)-sigm(1);
        inflpts(end+1) = indmed;
        steepnesses(end+1) = x(4);
        sigms(end+1,:) = sigm;
    end
    clf
    plot(data(ind,:))
    hold on
    plot(sigm,'g')
    hold off
    pause()
end

Pars.Ampls = ampls;
Pars.Baselines = baselines;
Pars.SigmAmpls = sigmampls;
Pars.Inflpts = inflpts;
Pars.Steepnessses = steepnesses;
Pars.Endlines = endlines;
Pars.Sigms = sigms;
