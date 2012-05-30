function TempSimData = SimulateSigmoidDetermEpid(Data,Parameters,HIVModel)

    tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
    Parameters.NbTSteps = length(tis);


    Parameters = HIV_Initialize(Parameters);
    
    rate = Parameters.Sigmrate.Value;
    base = Parameters.Sigmbase.Value;
    mu = Parameters.Sigmmu.Value;
    tinfl = Parameters.Sigmtinfl.Value;

    c = 1/(1+exp(tinfl/rate));
    b = (mu-base)*c/(1-c);
    a = base - b;
    
    
%     tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
%     Parameters.NbTSteps = length(tis);
%     Fts = baseline + Sigmoid((tis-tinfl)/rate)*(mu-baseline);
    
    tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
    Parameters.NbTSteps = length(tis);
    xis = [];
    xis(1) = exp(tinfl/rate);
    for i = 2:length(tis)
        xis(i) = xis(i-1) -1/rate*xis(i-1)*Parameters.ComputationTStep ;
    end
    Fts =  a + b./(c*(1+xis));
    
        

    
    
%     if not(isreal(Fts))
%         'r'
%     end
    
%     if or(sum(isnan(Fts)),sum(not(isreal(Fts))))
%         Fts = Fts(1)*ones(size(Fts));
%     end        
    
    Parameters.MultNoise = 0;
    TempSimData = HIV_CreateData(Fts,Parameters,HIVModel,Data);
    
    
    
    
    
    
    