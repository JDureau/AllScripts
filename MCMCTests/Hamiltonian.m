function H = Hamiltonian(x,p,Parameters)

M = Parameters.ScalingCov;
Dim = Parameters.Dim;

H = -log(max(Parameters.f(x,Parameters),eps))+1/2*log((2*pi)^Dim*det(M)) + 1/2*p'*M^(-1)*p;




