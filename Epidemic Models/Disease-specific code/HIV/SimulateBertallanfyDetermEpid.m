function TempSimData = SimulateBertallanfyDetermEpid(Data,Parameters,HIVModel)

    tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
    Parameters.NbTSteps = length(tis);


    Parameters = HIV_Initialize(Parameters);
    
    mu = Parameters.BRmu.Value;
    m = Parameters.BRmm1.Value + 1;
    baseline = Parameters.BRbase.Value;
    k = Parameters.k;
    
    
    tis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
    Parameters.NbTSteps = length(tis);
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
    
    
    
    
    
    
    