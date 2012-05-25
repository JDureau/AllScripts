function Gompertz = SimulateGompertzDetermEpid(Data,Parameters,HIVModel)

    Parameters = HIV_Initialize(Parameters);
    
    mu = Parameters.Gompmu.Value;
    baseline = Parameters.Gompbase.Value;
    inflpt = Parameters.Gomptinfl.Value;

    
    tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
    xis = [];
    xis(1) = (baseline^(1-m)-mu^(1-m))/(1-m);
    for i = 2:length(tis)
        xis(i) = xis(i-1) -k*xis(i-1)*Parameters.ComputationTStep ;
    end
    if baseline>m^(1/(1-m))*mu
        Fts =  ones(size(xis));
        Fts(2) = i;
    else
        Fts =  ((1-m)*xis+mu^(1-m)).^(1/(1-m));
    end
%     if not(isreal(Fts))
%         'r'
%     end
    
%     if or(sum(isnan(Fts)),sum(not(isreal(Fts))))
%         Fts = Fts(1)*ones(size(Fts));
%     end        
    
    Parameters.MultNoise = 0;
    TempSimData = HIV_CreateData(Fts,Parameters,HIVModel,Data);
    
    
    
    
    
    
    