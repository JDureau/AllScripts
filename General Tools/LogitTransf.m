function res = LogitTransf(Vect,Min,Max)

if nargin == 1
    Min = 0;
    Max = 1;
end

if not(and(Min-Vect<=0.0001*Min,Max-Vect>=-0.0001*Max)) 
    % Logit transf: out of bounds
    die
end

res = min(max(-10^100,log((Vect-Min)./(Max-Vect))),10^100);