function Pars = FitSigmoid(data,Parameters)

xinit(1) = data(1); %baseline
xinit(2) = max(data)-min(data); % ampl
[b,ind] = max(diff(data));
xinit(3) = ind; % inflpt
xinit(4) = 50;%b*Parameters.ComputationTStep;

[x,fval,exitflag,output] = fminsearch(@(x) DistToSigm(x,data,Parameters),xinit,optimset('MaxIter',1000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',500));

xis = (1:length(data))*Parameters.ComputationTStep;
sigm = x(1)+x(2)*Sigmoid((xis(1:length(data))-x(3)*Parameters.ComputationTStep)/x(4));


Pars.Baseline = sigm(1);
Pars.Endline = sigm(min(length(sigm),583));
Pars.SigmAmpl = x(2);
Pars.Ampl = sigm(end)-sigm(1);
Pars.Inflpt = x(3);
Pars.Steepness = x(4);
Pars.Sigm = sigm;
