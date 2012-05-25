function mLogLik = ComputemLogLik(PotValue,Parameters)

mLogLik = -eval(Parameters.LogLikFun);