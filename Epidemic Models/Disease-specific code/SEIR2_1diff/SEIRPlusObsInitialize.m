function Parameters = SEIRPlusObsInitialize(Parameters)

Parameters.EInitProp.Value = min(1,Parameters.EInitProp.Value);
Parameters.IInitProp.Value = min(1-Parameters.EInitProp.Value,Parameters.IInitProp.Value);
Parameters.RInitProp.Value = min(1-Parameters.EInitProp.Value-Parameters.IInitProp.Value,Parameters.RInitProp.Value);

Parameters.EInitProp.Value = max(0,Parameters.EInitProp.Value);
Parameters.IInitProp.Value = max(0,Parameters.IInitProp.Value);
Parameters.RInitProp.Value = max(0,Parameters.RInitProp.Value);

Parameters = UpdateParsNoTransfToTransf(Parameters);

InitialState(2) = Parameters.TotalPopulation*Parameters.EInitProp.Value;
InitialState(3) = Parameters.TotalPopulation*Parameters.IInitProp.Value;
InitialState(4) = Parameters.TotalPopulation*Parameters.RInitProp.Value;
InitialState(1) = max(0,Parameters.TotalPopulation -InitialState(2)-InitialState(3)-InitialState(4));
InitialState(5) = 0;
InitialState(6) = log(Parameters.betainit.Value);
InitialState(7) = log(Parameters.ObsPropinit.Value);

Parameters.InitialState = InitialState';