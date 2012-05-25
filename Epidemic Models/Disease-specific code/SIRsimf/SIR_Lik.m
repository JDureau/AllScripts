function Lik = SIR_Lik(Variables,Observation,Parameters)

Lik = normpdf(Variables(:,2),Observation,Observation*Parameters.SigmaObs);
