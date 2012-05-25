function Parameters = SIRsimple_Initialize(Parameters)

p_t = Parameters.p_t;

Parameters.InitialState(1) = Parameters.Sinit.Value*p_t;
Parameters.InitialState(3) = Parameters.Iinit.Value*p_t;
Parameters.InitialState(9) = log(Parameters.r0_init.Value);




