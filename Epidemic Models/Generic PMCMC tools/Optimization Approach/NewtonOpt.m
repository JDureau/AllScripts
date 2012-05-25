function res = NewtonOpt(Function,Parameters,Data,x0)

NbParticules = 15000;
it = 1;
x = x0;
[f,g,h] = Function(x, Parameters, Data, NbParticules);

xvals = x0;
fvals = f;
gvals = g;
hvals = h;
steps = -g/h;
while it <100
    it = it+1;
    disp(it)
    
    x = x - g/h;
    [f,g] = Function(x, Parameters, Data, NbParticules);
    xvals(end+1) = x;
    fvals(end+1) = f;
    gvals(end+1) = g;
    hvals(end+1) = h;
    steps(end+1) = -g/h;
end

res.x = x;
res.fval = f;
res.grad = g;
    