function Lambda = ComputeLambda(nbpoints_discr,step,H)

% create first row of C
tmp = 0:(nbpoints_discr-1);
v = step^(-2*H)*(0.5*(abs(tmp+1)).^(2*H) + 0.5*(abs(tmp-1)).^(2*H) - (abs(tmp)).^(2*H));
v(1) = step^(-2*H);
Lambda = diag(fft([v v(end) v(end:-1:2)]));