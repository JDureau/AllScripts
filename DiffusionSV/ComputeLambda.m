function Lambda = ComputeLambda(N,step,H)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;


% create first row of C
tmp = 0:(N-2);
v = step^(-2*H)*(0.5*(abs(tmp+1)).^(2*H) + 0.5*(abs(tmp-1)).^(2*H) - (abs(tmp)).^(2*H));
v(1) = step^(-2*H);
Lambda = diag(fft([v v(end) v(end:-1:2)]));