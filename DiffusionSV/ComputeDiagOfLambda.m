function DiagOfLambda = ComputeDiagOfLambda(N,step,H)
    % Length of Z : 2*(N)
    % Length of Bh : N
    % Length of X : N  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


% create first row of C
tmp = 0:(N-1);
v = step^(2*H)*(0.5*(abs(tmp+1)).^(2*H) + 0.5*(abs(tmp-1)).^(2*H) - (abs(tmp)).^(2*H));
v(1) = step^(2*H);
% DiagOfLambda = (fft([v v(end) v(end:-1:2)]));

DiagOfLambda = (fft([v 0 v(end:-1:2)]));

