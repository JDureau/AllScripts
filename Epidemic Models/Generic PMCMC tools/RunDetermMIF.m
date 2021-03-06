function Res = RunDetermMIF(Data, Model, Parameters)


L = Parameters.MIFL;
NbIts = Parameters.MIFNbIts;
a = Parameters.MIFa;
b = Parameters.MIFb;


Parameters.RunningMif = 1;
Names = Parameters.Names.Estimated;
NbCstPars = length(Names);
NbInit = length(Parameters.Names.EstimatedInit);
NbNotInit = NbCstPars - NbInit;
Parameters.MIFNbInit = NbInit;
Parameters.MIFNbNotInit = NbNotInit;

Parameters.MIFVarInit = 1;
Parameters.MIFVarNotInit = 1;

Thetas = [];
TransfThetas = [];
LogLiks = [];
for IndIt = 1:NbIts
    Parameters.MIFVarInit    = a^(IndIt - 1)*b*Parameters.MIFVarInit;
    Parameters.MIFVarNotInit = a^(IndIt - 1)*Parameters.MIFVarNotInit;
    
    TmpRes = EstimationEKFGen(Data, Model, Parameters);
    LogLiks(IndIt) = TmpRes.LogLik;
    
    nbstatevars = size(TmpRes.PosteriorMeans,1) - NbCstPars;
    for  i = 1:NbCstPars
        if Parameters.(Names{i}).Init
            Parameters.(Names{i}).TransfValue = TmpRes.PosteriorMeans(nbstatevars +i,L);
        else
            Parameters.(Names{i}).TransfValue = TmpRes.PosteriorCovs(1,nbstatevars +i,nbstatevars +i)*(diff(TmpRes.PosteriorMeans(nbstatevars +i,:)))*((squeeze(TmpRes.PosteriorCovs(2:end,nbstatevars +i,nbstatevars +i))).^(-1));
        end
    end
    Parameters = UpdateParsTransfToNoTransf(Parameters);
    for  i = 1:NbCstPars
        Thetas(i,IndIt) = Parameters.(Names{i}).Value;
        TransfThetas(i,IndIt) = Parameters.(Names{i}).TransfValue;
    end
    disp(['IndIt = ' num2str(IndIt) ', LogLik = ' num2str(TmpRes.LogLik)])

    tmpThetas = TmpRes.PosteriorMeans(nbstatevars+1:end,:);
    k = ceil(sqrt(length(Parameters.Names.Estimated)));
    for  i = 1:NbCstPars
        subplot(k,k,i)
        plot(tmpThetas(i,:));
        title(Names{i});
    end
end

Res.Thetas = Thetas;
Res.TransfThetas = TransfThetas;
Res.Parameters = Parameters;
Res.LogLiks = LogLiks;


    