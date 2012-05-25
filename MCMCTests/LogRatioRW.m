function LogRatio = LogRatioRW(x,xstar,Parameters)


% numerator (xstar)
LogLik = max(log(eps),log(Parameters.f(xstar',Parameters)));
LogNum = LogLik;

% denominator (x)
LogLik = max(log(eps),log(Parameters.f(x',Parameters)));
LogDenom = LogLik;

LogRatio = LogNum - LogDenom;

if isnan(LogRatio)
    disp('LogRatio problem')
end

