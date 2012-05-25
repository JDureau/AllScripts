function Parameters = SIRsimf_Initialize(Parameters)

p_t = Parameters.p_t;

Parameters.InitialState(1) = Parameters.Sinit.Value*p_t;
Parameters.InitialState(2) = Parameters.Sinit.Value*p_t;
Parameters.InitialState(3) = Parameters.Iinit.Value*p_t/2;
Parameters.InitialState(4) = Parameters.Iinit.Value*p_t/2;
Parameters.InitialState(5) = Parameters.Iinit.Value*p_t/2;
Parameters.InitialState(6) = Parameters.Iinit.Value*p_t/2;
Parameters.InitialState(7) = 0;
Parameters.InitialState(8) = 0;
Parameters.InitialState(9) = log(Parameters.r0_init.Value);




