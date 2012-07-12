function LogRatio = LogRatioGMHMC(bigx,bigxstar,Parameters)

Epsil = Parameters.Epsil;
epsilon = 10^(-10);
Dim = Parameters.Dim;

x = bigx(1:Dim)';
xstar = bigx(1:Dim)';
pstar = bigxstar(Dim+1:2*Dim)';
pnp1 = bigxstar(2*Dim+1:3*Dim)';

Hstar = Hamiltonian(xstar,pstar,Parameters);
H = Hamiltonian(x,pnp1,Parameters);
LogRatio = -(H - Hstar);
