function Pars = SampleUnif( ParsInf, ParsSup)

Pars = rand(1,length(ParsInf))*(ParsSup-ParsInf) + ParsInf;