function res = logit(x,Min,Max)

if nargin == 1
    Min = 0;
    Max = 1;
end

if not(and(Min-x<=0.0001*Min,Max-x>=-0.0001*Max)) 
    % Logit transf: out of bounds
    die
end

res = min(max(-10^100,log((x-Min)./(Max-x))),10^100);