function Parameters = SIRsimple_Initialize(Parameters)

p_t = Parameters.p_t;

Parameters.InitialState(1) = Parameters.Sinit.Value*p_t;
Parameters.InitialState(2) = Parameters.Iinit.Value*p_t;
Parameters.InitialState(3) = 0;
Parameters.InitialState(4) = log(Parameters.r0_init.Value);




