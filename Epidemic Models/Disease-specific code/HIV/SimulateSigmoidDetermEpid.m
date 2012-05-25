function TempSimData = SimulateSigmoidDetermEpid(Data,Parameters,HIVModel)

    tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
    Parameters.NbTSteps = length(tis);


    Parameters = HIV_Initialize(Parameters);
    
    mu = Parameters.Sigmmu.Value;
    rate = Parameters.Sigmrate.Value;
    baseline = Parameters.Sigmbase.Value;
    tinfl = Parameters.Sigmtinfl.Value;
    
    
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
    Fts =  baseline + 1./(1+xis)*(mu-baseline);
    
        

    
    
%     if not(isreal(Fts))
%         'r'
%     end
    
%     if or(sum(isnan(Fts)),sum(not(isreal(Fts))))
%         Fts = Fts(1)*ones(size(Fts));
%     end        
    
    Parameters.MultNoise = 0;
    TempSimData = HIV_CreateData(Fts,Parameters,HIVModel,Data);
    
    
    
    
    
    
    