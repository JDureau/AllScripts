function Ress = Burnin(Parameters,N)

Parameters.CholCov = chol(2.38^2/Parameters.Dim*(-Parameters.Hess)^-1);
Parameters.aim = 0.23;
Dim = Parameters.Dim;

Ress = {};

InitValue = Parameters.ArgMax;
Ress{1} = RunMCMC(InitValue,Parameters,10000);
plot(Ress{1}.Vals(1,:),Ress{1}.Vals(2,:),'.')

Parameters.aim = 0.6;

for i = 2:N
    Parameters.aim = 0.6;
    if strcmp(Parameters.TypeMeth,'RW')
        Parameters.CholCov = chol(cov(Ress{i-1}.Vals'));
    elseif strcmp(Parameters.TypeMeth,'GMM')
        tmps = {};
        BICs = [];
        AICs = [];
        clf
        for j = 1:10
            tmps{j} = gmdistribution.fit(Ress{i-1}.Vals',j);
            lnL = sum(log(pdf(tmps{j},Ress{i-1}.Vals')));
            k = Dim*j+j-1+j*Dim*(Dim-1)/2;
            BICs(j) = -2*lnL + k*log(size(Ress{i-1}.Vals,2));
            scatter(Ress{i-1}.Vals(1,:),Ress{i-1}.Vals(2,:),'.')
            hold on
            MinX = min(Ress{i-1}.Vals(1,:));
            MaxX = max(Ress{i-1}.Vals(1,:));
            MinY = min(Ress{i-1}.Vals(2,:));
            MaxY = max(Ress{i-1}.Vals(2,:));
            
            
            ezcontour(@(x,y)pdf(tmps{j},[x y]),[MinX MaxX],[MinY MaxY]);
            hold off
            title([BICs(j)])
            pause(0.01)
        end
        [b,ind] = min(BICs);
        Parameters.Dens = tmps{ind};
    end
    [Parameters, TempPar] = CalibrateMCMC( Parameters.ArgMax , Parameters); 
    Ress{i} = RunMCMC(InitValue,Parameters,10000);
    plot(Ress{i}.Vals(1,:),Ress{i}.Vals(2,:),'.')
    pause(0.01)
end
