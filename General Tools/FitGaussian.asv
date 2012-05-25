function Fit = FitGaussian(xvect,yvect,Parameters)
ampl = max(xvect') - min(xvect');
MinVect = min(xvect') - 0.1*ampl;
MaxVect = max(xvect') + 0.1*ampl;

Initialization = [];
Initialization(1) = mean(exp(yvect))^-(1);
c = Initialization(1);
MinC = c/10000000000000000000000000000;
MaxC = c*10000000000000000000000000000;
Initialization(1) = LogitTransf(c,MinC,MaxC);
dim = size(xvect,1);
% [bof,ind] = max(yvect);
% Initialization(2:1+dim) = LogitTransf(xvect(:,ind)',MinVect,MaxVect);
for i = 1:dim
    Initialization(1+i) = LogitTransf((xvect(i,:)*exp(yvect)')/sum(exp(yvect))',MinVect(i),MaxVect(i));
end

% temp = triu(10*max(ampl.^2)*ones(dim,dim)) + diag((ampl*2).^2);
% MaxBounds = reshape(nonzeros(temp),dim*(dim+1)/2,1);
% MinBounds = zeros(size(MaxBounds));
temp = reshape(nonzeros(chol(cov(xvect'))),dim*(dim+1)/2,1);
Initialization(dim+2:dim+1+dim*(dim+1)/2) = temp;
inds = ceil(rand(1,min(200,ceil(size(xvect,2)*0.1)))*size(xvect,2));
[x,fval,exitflag,output] = fminsearch(@(x)TestNormalModel(x,xvect',yvect,MinVect,MaxVect,MinC,MaxC,ampl,inds),Initialization,optimset('MaxFunEval',10^6,'MaxIter',10^6,'TolX',1e-8,'TolFun',1e-8));
a = x(dim+2:end); 
b = triu(ones(dim),1)+eye(dim); 
b(~~b)=a ;
Cov = b*b';
c = InvLogitTransf(x(1),MinC,MaxC);
mu = ones(1,dim);
for i = 1:dim
    mu(i) = InvLogitTransf(x(1+i),MinVect(i),MaxVect(i));
end
disp(fval)
disp(output.message)
disp(mean(exp(yvect))^-(1)/c)

figure(1)
clf
k = ceil(sqrt(dim));
xvalues = xvect';
yvalues = yvect;
for i = 1:dim
    subplot(k,k,i) 
    plot(xvalues(:,i),log(c*exp(yvalues)),'.')
    hold on
    perc = 5;
    inds = ceil(rand(1,floor(length(xvalues(:,i))/perc))*length(xvalues(:,i)));
    yis = mvnpdf(xvalues(inds,:),mu,Cov);
    plot(xvalues(inds,i),log(yis),'+g')
    hold off
end
figure(2)
clf
k = ceil(sqrt(dim));
for i = 1:dim
    subplot(k,k,i) 
    plot(xvalues(:,i))
    hold on
    xis = 1:length(xvalues(:,i));
    plot(xis,mu(i)*ones(size(xis)),'k')
    plot(xis,(mu(i)+sqrt(Cov(i,i)))*ones(size(xis)),'r')
    plot(xis,(mu(i)-sqrt(Cov(i,i)))*ones(size(xis)),'r')
    try
        plot(xis,Parameters.RealDens.mu(i)*ones(size(xis)),'g')
        plot(xis,(Parameters.RealDens.mu(i)+sqrt(Parameters.RealDens.Sigma(i,i)))*ones(size(xis)),'--g')
        plot(xis,(Parameters.RealDens.mu(i)-sqrt(Parameters.RealDens.Sigma(i,i)))*ones(size(xis)),'--g')
    end
    hold off
end


Fit.Cov = Cov;
Fit.Mu = x(2:dim+1);
Fit.c = x(1);

