function xstar = SampleRW(x,Parameters)

xstar = x+Parameters.Epsil*Parameters.CholCov*randn(Parameters.Dim,1);