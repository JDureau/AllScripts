function Parameters = FindFisherInfMat(CurrentPos,Parameters) 

Test = 0;
n = 10000;
NbItsMax = 50;
NbIts = 0;

while not(Test)
    [x,fval,exitflag,output] = fminsearch(@(x) -Parameters.f(x,Parameters),CurrentPos,optimset('MaxIter',n,'TolX',1e-8,'TolFun',1e-7));
    CurrentPos = x;
    Parameters = ComputeHess(CurrentPos,Parameters);
    if mean(eig(-Parameters.Hess)>0)==1
        Test = 1;
    end
    NbIts = NbIts + 1;
    if NbIts >NbItsMax
        Test = 1;
        disp('Could''nt maximize')
    end
end
CurrentPos
Parameters.ArgMax = CurrentPos;
disp(Parameters.Hess);
