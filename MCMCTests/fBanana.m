function f = fBanana(x,Parameters)
% from Robers and Rosenthal 08
n = length(x);
B = Parameters.B;

tmp = -x(1)^2/200 - 0.5*(x(2)+B*x(1)^2-100*B)^2;
if length(x)>2
    tmp = tmp -0.5*sum((x(3:end)).^2);
end

f = exp(tmp);
if or(isnan(f),isinf(f))
    disp('Banana problem')
end


