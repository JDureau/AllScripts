function Result = OptParEstimation(FilteringFunction, Data, Parameters, x0)

mLogLikInit = FilteringFunction(x0, Data, Parameters, Parameters.NbParticules, 'ParEst');
history.x = [];
history.fval = [];
searchdir = [];


A = [1 0 ; 0 1; -1 0; 0 -1];
b = [Parameters.BetaSup; Parameters.GammaSup; -Parameters.BetaInf; -Parameters.GammaInf];
options = optimset('outputfcn',@outfun,'display','off');
warning off
[OptPars,OptmLogLik,flag] = fmincon(@(x) FilteringFunction(x, Data, Parameters, Parameters.NbParticules, 'ParEst'),x0,A,b,[],[],[],[],[],options);

Result.InitPars = x0;
Result.InitLogLik = exp(-mLogLikInit);
Result.InitLik = exp(-mLogLikInit);
Result.OptPars = OptPars;
Result.OptLogLik = -OptmLogLik;
Result.OptLik = exp(-OptmLogLik);
Result.Flag = flag;
Result.HistoryX = history.x;
Result.HistoryFval = -history.fval;

 
 function stop = outfun(x,optimValues,state)
     stop = false;
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
         % Concatenate current search direction with 
%          % searchdir.
%            searchdir = [searchdir;... 
%                         optimValues.searchdirection'];
%            plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
%            text(x(1)+.15,x(2),... 
%                 num2str(optimValues.iteration));
%            title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
 end


end