function Parameters = HIVECP_Initialize(Parameters)

InitialState = [];

InitialState(1) = 1-Parameters.InitialI.Value;
InitialState(2) = Parameters.InitialI.Value;
InitialState(3) = 0;
InitialState(4) = Parameters.InitialI.Value;

Parameters.InitialState = InitialState';

