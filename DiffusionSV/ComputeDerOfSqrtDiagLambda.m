function DerOfSqrtDiagLambda = ComputeDerOfSqrtDiagLambda(N,step,H)
    % Length of Z : 2*(N-1)
    % Length of Bh : N-1
    % Length of X : N-1  
    %    X(i) = X(i-1) + Bh(i) with X(0) = 0;

    
% create first row of C

tmp = 0:(N-1);
v = step^(-2*H)*(0.5*(abs(tmp+1)).^(2*H) + 0.5*(abs(tmp-1)).^(2*H) - (abs(tmp)).^(2*H));
v(1) = step^(-2*H);
DiagOfLambda = (fft([v v(end) v(end:-1:2)]));


tmp = 0:(N-1);

v = 2*log(step)*step^(2*H)*(0.5*(abs(tmp+1)).^(2*H) + 0.5*(abs(tmp-1)).^(2*H) - (abs(tmp)).^(2*H)) + ...
     step^(2*H)*(log(abs(tmp+1)).*(abs(tmp+1)).^(2*H) + log(abs(tmp-1)).*(abs(tmp-1)).^(2*H) - 2*log(abs(tmp)).*(abs(tmp)).^(2*H));

v(1) = 2*log(step)*step^(2*H);

v(2) = 2*log(step)*step^(2*H)*(0.5*(2).^(2*H) - 1) + step^(2*H)*(log(2).*(2).^(2*H));



DiagOfDerLambda = (fft([v v(end) v(end:-1:2)]));


DiagOfLambda = ComputeDiagOfLambda(N,step,H);


DerOfSqrtDiagLambda = (DiagOfDerLambda./(2*sqrt((DiagOfLambda))));


