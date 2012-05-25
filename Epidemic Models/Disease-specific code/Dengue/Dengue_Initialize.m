function Parameters = Dengue_Initialize(Parameters)


InitialState(2) = Parameters.I1Init.Value*Parameters.PopulationSize.Value;
InitialState(3) = Parameters.I2Init.Value*Parameters.PopulationSize.Value;
InitialState(4) = Parameters.Q1Init.Value*Parameters.PopulationSize.Value;
InitialState(5) = Parameters.Q2Init.Value*Parameters.PopulationSize.Value;
InitialState(6) = Parameters.R1Init.Value*Parameters.PopulationSize.Value;
InitialState(7) = Parameters.R2Init.Value*Parameters.PopulationSize.Value;
InitialState(8) = Parameters.I12Init.Value*Parameters.PopulationSize.Value;
InitialState(9) = Parameters.I21Init.Value*Parameters.PopulationSize.Value;
InitialState(1) = Parameters.SInit.Value*Parameters.PopulationSize.Value;
InitialState(10) = Parameters.r0Init.Value;
InitialState(11) = 0;

Parameters.InitialState = InitialState';