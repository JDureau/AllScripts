function xopt = FindOptimalValue(InitValue,Parameters)

[xopt,fval,exitflag,output] = fminsearch(@(x) ComputemLogLik(x,Parameters),InitValue,optimset('MaxIter',4000,'TolX',1e-8,'TolFun',1e-8,'MaxFunEvals',10000));

disp(output.message)