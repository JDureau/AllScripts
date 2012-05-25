function mpred = SEIRpredict(mpred,Parameters)

    mpred(1:5,1) = max(mpred(1:5,1),0);
    mpred(1:5,1) = min(mpred(1:5,1),Parameters.TotalPopulation);
    mtemp = mpred;
    
    beta = exp(mtemp(6));
    mpred(1) = mpred(1) + (-beta*mtemp(1)*mtemp(3)/TotPop)*TStep;
    mpred(2) = mpred(2) + ( beta*mtemp(1)*mtemp(3)/TotPop- Parameters.k.Value*mtemp(2))*TStep;
    mpred(3) = mpred(3) + ( Parameters.k.Value*mtemp(2) - Parameters.gamma.Value*mtemp(3))*TStep;
    mpred(4) = mpred(4) + ( Parameters.gamma.Value*mtemp(3))*TStep;
    mpred(5) = mpred(5) + ( Parameters.k.Value*mtemp(2))*TStep;
    mpred(6) = mtemp(6);