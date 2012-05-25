function ResGen = GenerateDataForGivenSig(Data,Parameters,Model,Fts)



xis = Parameters.ComputationTStep:Parameters.ComputationTStep:Parameters.ObservationLength;
Parameters.NbTSteps = length(xis);

% Names = Parameters.Names.All;
% for i = 1:length(Names)
%     Parameters.(Names{i}).Sample = 1;
% end
it = 0;
test = 0;
while not(test)
    ampl = rand(1,1);
    baseline = rand(1,1)*(1-ampl);
    asympt = baseline+ampl;
    
    
    Parameters = SampleParameters(Parameters);

    inflpt = ceil(min(300,max(0,abs(randn(1,1)*150))));
    steepness = rand(1,1)*50;

    if nargin == 3
        Fts = baseline+ampl*Sigmoid((xis-inflpt)/steepness);
    end
        
    Parameters.InitialFt.Value = Fts(1);    
    Parameters = UpdateParsNoTransfToTransf(Parameters);
    
    DataGen = HIV_CreateData(Fts,Parameters,Model,Data);
    if and(DataGen.BuiltTraj(600,7)>2,DataGen.BuiltTraj(600,7)<40)
        if and(DataGen.BuiltTraj(600,8)>2,DataGen.BuiltTraj(600,8)<40)
            test = 1;
        end
    end
%     test 
    try
        Parameters.InitialFt.Value = DataGen.BuiltTraj(1,9);
        Parameters.InitialFt.Estimated;
        Parameters = DefineEstimatedParametersIndexes(Parameters);
        Parameters = DefineTransfFunctions(Parameters);
        Parameters = DefinePriors(Parameters);
        Parameters = UpdateParsNoTransfToTransf(Parameters);
    catch
        test = 0;
    end
    it = it+1;
    if it>25
        test = 1;
        disp('crash')
        die
    end
%     if rand(1)>DataGen.BuiltTraj(end,9)-DataGen.BuiltTraj(1,9)
%         test = 0;
%     end
end
DataGen.Fts = Fts;
ResGen.ampl = ampl;
ResGen.infl = inflpt;
ResGen.Data = DataGen;
ResGen.Parameters = Parameters;