function f = fGMM(x,Parameters)

if and(size(x,1) == 1,size(x,2) == Parameters.Dim)
    f = pdf(Parameters.RealDens,x);
elseif and(size(x,2) == 1,size(x,1) == Parameters.Dim)
    f = pdf(Parameters.RealDens,x');
else
    'Wrong argument in fGMM'
    die
end
