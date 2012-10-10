function Parameters = SEIRInitialize2(Parameters)

die
InitialState(1) = log(Parameters.TotalPopulation);
InitialState(2) = log(Parameters.TotalPopulation*Parameters.EInitProp.Value);
InitialState(3) = log(Parameters.TotalPopulation*Parameters.IInitProp.Value);
InitialState(4) = log(Parameters.TotalPopulation*Parameters.RInitProp.Value);
InitialState(5) = log(eps);
InitialState(6) = log(Parameters.betainit.Value);

Parameters.InitialState = InitialState';