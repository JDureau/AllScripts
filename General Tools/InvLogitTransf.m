function res = InvLogitTransf(Vect,Min,Max)


if nargin == 1
    Min = 0*ones(size(Vect));
    Max = 1*ones(size(Vect));
end


res = (Max.*exp(Vect)+Min)./(1+exp(Vect));
inds = find(or(isnan(res),isinf(res)));
res(inds) = (Max*ones(size(inds))+Min*ones(size(inds)).*exp(-Vect(inds)))./(1+exp(-Vect(inds)));
    
if not(and(res>=Min-0.001*abs(Min),res<=Max+0.001*abs(Max)))
    % Logit transf: out of bounds
    die
end
