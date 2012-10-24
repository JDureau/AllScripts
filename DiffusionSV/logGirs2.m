function [log_gir] = logGirs2(Q,MU,step)

%Returns the log Girsanov density for a GCIR diffusion transformed to Xdot
%MU is the vector with the drift values

%stochastic integral
stoc_int1 = MU(:,1:end-1).*diff(Q,1,2); 
stoc_int = sum(sum(stoc_int1));

%path integral
int = sum(sum(MU.^2))*step;

log_gir = (stoc_int - 0.5*int);


