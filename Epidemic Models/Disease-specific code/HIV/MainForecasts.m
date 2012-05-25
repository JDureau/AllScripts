%% load paths

cd('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Filtering')
addpath('H:\My Documents\PhD Work\Matlab Scripts\General Tools')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes\Resampling\pf_resampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Joint Sampling')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\MIF')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Model Selection')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Optimization Approach')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\Models')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Toolboxes')
addpath('H:\My Documents\PhD Work\Matlab Scripts\Epidemic Models\Generic Parameter Estimation\HIV')


% Main HIV Forecasts

SavePath = 'S:\Results';
Res = load([SavePath '\HIV_Shimoga.mat' ]);
Res  = Res.ResRW;
Parameters = Res.Parameters;
Names = Parameters.Names.Estimated;

NbYearsForecast = 10;
NbTimeSteps = NbYearsForecast*12/Parameters.ComputationTStep;
NbIts = size(Res.Paths,1);

record = zeros(NbIts,9,NbTimeSteps);

go_back = 0;
ParsVect = Res.Thetas;
mpred = squeeze(Res.Paths(:,:,size(Res.Paths,3)-go_back));
Fts = squeeze(Res.Paths(:,9,size(Res.Paths,3)-go_back));
record = [];
TStep = Parameters.ComputationTStep;
BetasFM = ones(1,size(mpred,1))-(1-ParsVect(Parameters.BetaFMPerAct.Index,:)).^ParsVect(Parameters.NumberActsPerClient.Index,:);
BetasMF = ones(1,size(mpred,1))-(1-ParsVect(Parameters.BetaMFPerAct.Index,:)).^ParsVect(Parameters.NumberActsPerClient.Index,:);
for IndDiscr = 1:NbTimeSteps;
    mtemp = mpred;
    TotF1 = mtemp(:,1) + mtemp(:,2);
    TotF2 = mtemp(:,3) + mtemp(:,4);
    TotM  = mtemp(:,5) + mtemp(:,6);
    mtemp(:,9) = Fts;
    mpred(:,1) = mpred(:,1) + ( (ParsVect(Parameters.MuFm1.Index,:).^-1)'.*(TotF1) + (ParsVect(Parameters.Alpham1.Index,:).^-1)'.*mtemp(:,2) - BetasMF'.*ParsVect(Parameters.CF1.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,1).*mtemp(:,6)./TotM-(ParsVect(Parameters.MuFm1.Index,:).^-1)'.*mtemp(:,1))*TStep;
    mpred(:,2) = mpred(:,2) + ( BetasMF'.*ParsVect(Parameters.CF1.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,1).*mtemp(:,6)./TotM - ((ParsVect(Parameters.MuFm1.Index,:).^-1)' + (ParsVect(Parameters.Alpham1.Index,:).^-1)').*mtemp(:,2))*TStep;
    mpred(:,3) = mpred(:,3) + ( (ParsVect(Parameters.MuFm1.Index,:).^-1)'.*(TotF2) + (ParsVect(Parameters.Alpham1.Index,:).^-1)'.*mtemp(:,4) - BetasMF'.*ParsVect(Parameters.CF2.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,3).*mtemp(:,6)./TotM-(ParsVect(Parameters.MuFm1.Index,:).^-1)'.*mtemp(:,3))*TStep;
    mpred(:,4) = mpred(:,4) + ( BetasMF'.*ParsVect(Parameters.CF2.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,3).*mtemp(:,6)./TotM - ((ParsVect(Parameters.MuFm1.Index,:).^-1)' + (ParsVect(Parameters.Alpham1.Index,:).^-1)').*mtemp(:,4))*TStep;
    mpred(:,5) = mpred(:,5) + ( (ParsVect(Parameters.MuMm1.Index,:).^-1)'.*(TotM) + (ParsVect(Parameters.Alpham1.Index,:).^-1)'.*mtemp(:,6) - BetasFM'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,5).*(ParsVect(Parameters.CF1.Index,:)'.*TotF1./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,2)./TotF1+ParsVect(Parameters.CF2.Index,:)'.*TotF2./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,4)./TotF2)-(ParsVect(Parameters.MuMm1.Index,:).^-1)'.*mtemp(:,5))*TStep;
    mpred(:,6) = mpred(:,6) + ( BetasFM'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,5).*(ParsVect(Parameters.CF1.Index,:)'.*TotF1./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,2)./TotF1+ParsVect(Parameters.CF2.Index,:)'.*TotF2./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,4)./TotF2) - ((ParsVect(Parameters.MuMm1.Index,:).^-1)'+(ParsVect(Parameters.Alpham1.Index,:).^-1)').*mtemp(:,6))*TStep;
    mpred(:,7) = (mpred(:,2) + mpred(:,4))./(TotF1+TotF2)*100; 
    mpred(:,8) = (mpred(:,6))./(TotM)*100; 
    record(:,:,IndDiscr) = mpred;
end
    


PlotForecast(Res,record,go_back)

% Loop

Names = {'\HIV_Mysore_3rounds','\HIV_Belgaum_3rounds','\HIV_Bellary','\HIV_EastGodavry','\HIV_Guntur','\HIV_Hyderabad','\HIV_Yevatmal','\HIV_Shimoga'};

for i = 1:length(Names)
    SavePath = 'S:\Results';
    Res = load([SavePath Names{i} ]);
    Res = Res.ResRW;
    Parameters = Res.Parameters;


    NbYearsForecast = 41;
    NbTimeSteps = NbYearsForecast*12/Parameters.ComputationTStep;
    NbIts = size(Res.Paths,1);

    
    year = floor(size(Res.Paths,3)/24) + 1984;
    
    clear  FSW_Prev M_Prev FSW_Cases M_Cases
    FSW_Prev{1} = {'Year','Mean','q(0.25)','q(0.755)','q(0.025)','q(0.975)'};
    M_Prev{1} = {'Year','Mean','q(0.25)','q(0.755)','q(0.025)','q(0.975)'};
    FSW_Cases{1} = {'Year','Mean','q(0.25)','q(0.755)','q(0.025)','q(0.975)'};
    M_Cases{1} = {'Year','Mean','q(0.25)','q(0.755)','q(0.025)','q(0.975)'};
    years = [];

    go_back = 0;
    ParsVect = Res.Thetas;
    mpred = squeeze(Res.Paths(:,:,size(Res.Paths,3)-go_back));
    Fts = squeeze(Res.Paths(:,9,size(Res.Paths,3)-go_back));
    record = zeros(NbIts,9,NbTimeSteps);
    TStep = Parameters.ComputationTStep;
    BetasFM = 1-(1-ParsVect(Parameters.BetaFMPerAct.Index,:)).^ParsVect(Parameters.NumberActsPerClient.Index,:);
    BetasMF = 1-(1-ParsVect(Parameters.BetaMFPerAct.Index,:)).^ParsVect(Parameters.NumberActsPerClient.Index,:);
    
    cpt = 0;
    NewCases_FSW = 0;
    NewCases_M = 0;
    for IndDiscr = 1:NbTimeSteps;
        mtemp = mpred;
        TotF1 = mtemp(:,1) + mtemp(:,2);
        TotF2 = mtemp(:,3) + mtemp(:,4);
        TotM  = mtemp(:,5) + mtemp(:,6);
        mtemp(:,9) = Fts;
        mpred(:,1) = mpred(:,1) + ( (ParsVect(Parameters.MuFm1.Index,:).^-1)'.*(TotF1) + (ParsVect(Parameters.Alpham1.Index,:).^-1)'.*mtemp(:,2) - BetasMF'.*ParsVect(Parameters.CF1.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,1).*mtemp(:,6)./TotM-(ParsVect(Parameters.MuFm1.Index,:).^-1)'.*mtemp(:,1))*TStep;
        mpred(:,2) = mpred(:,2) + ( BetasMF'.*ParsVect(Parameters.CF1.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,1).*mtemp(:,6)./TotM - ((ParsVect(Parameters.MuFm1.Index,:).^-1)' + (ParsVect(Parameters.Alpham1.Index,:).^-1)').*mtemp(:,2))*TStep;
        mpred(:,3) = mpred(:,3) + ( (ParsVect(Parameters.MuFm1.Index,:).^-1)'.*(TotF2) + (ParsVect(Parameters.Alpham1.Index,:).^-1)'.*mtemp(:,4) - BetasMF'.*ParsVect(Parameters.CF2.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,3).*mtemp(:,6)./TotM-(ParsVect(Parameters.MuFm1.Index,:).^-1)'.*mtemp(:,3))*TStep;
        mpred(:,4) = mpred(:,4) + ( BetasMF'.*ParsVect(Parameters.CF2.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,3).*mtemp(:,6)./TotM - ((ParsVect(Parameters.MuFm1.Index,:).^-1)' + (ParsVect(Parameters.Alpham1.Index,:).^-1)').*mtemp(:,4))*TStep;
        mpred(:,5) = mpred(:,5) + ( (ParsVect(Parameters.MuMm1.Index,:).^-1)'.*(TotM) + (ParsVect(Parameters.Alpham1.Index,:).^-1)'.*mtemp(:,6) - BetasFM'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,5).*(ParsVect(Parameters.CF1.Index,:)'.*TotF1./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,2)./TotF1+ParsVect(Parameters.CF2.Index,:)'.*TotF2./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,4)./TotF2)-(ParsVect(Parameters.MuMm1.Index,:).^-1)'.*mtemp(:,5))*TStep;
        mpred(:,6) = mpred(:,6) + ( BetasFM'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,5).*(ParsVect(Parameters.CF1.Index,:)'.*TotF1./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,2)./TotF1+ParsVect(Parameters.CF2.Index,:)'.*TotF2./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,4)./TotF2) - ((ParsVect(Parameters.MuMm1.Index,:).^-1)'+(ParsVect(Parameters.Alpham1.Index,:).^-1)').*mtemp(:,6))*TStep;
        mpred(:,7) = (mpred(:,2) + mpred(:,4))./(TotF1+TotF2)*100; 
        mpred(:,8) = (mpred(:,6))./(TotM)*100; 
        NewCases_FSW = NewCases_FSW + (BetasMF'.*ParsVect(Parameters.CF1.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,1).*mtemp(:,6)./TotM)*TStep + ( BetasMF'.*ParsVect(Parameters.CF2.Index,:)'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,3).*mtemp(:,6)./TotM)*TStep;
        NewCases_M = NewCases_M + (BetasFM'.*(1-ParsVect(Parameters.eHIV.Index,:)'.*mtemp(:,9)).*mtemp(:,5).*(ParsVect(Parameters.CF1.Index,:)'.*TotF1./(ParsVect(Parameters.CF2.Index,:)'.*TotF2+ParsVect(Parameters.CF1.Index,:)'.*TotF1).*mtemp(:,2)./TotF1))*TStep;
        record(:,:,IndDiscr) = mpred;
        cpt = cpt+1;
        if cpt == 24
            cpt = 0;
            year = year + 1;
            years(end+1) = year;
            
            ind = size(FSW_Prev,1)+1;
            
            tmp = mpred(:,7);
            tmp = {year mean(tmp) quantile(tmp,0.25) quantile(tmp,0.75) quantile(tmp,0.025) quantile(tmp,0.975)};
            for j = 1:6                
                FSW_Prev{ind,j} = tmp{j};
            end
            
            tmp = mpred(:,8);
            tmp = {year mean(tmp) quantile(tmp,0.25) quantile(tmp,0.75) quantile(tmp,0.025) quantile(tmp,0.975)};
            for j = 1:6                
                M_Prev{ind,j} = tmp{j};
            end
            
            tmp = NewCases_FSW./(TotF1+TotF2)*1000;
            tmp = {year mean(tmp) quantile(tmp,0.25) quantile(tmp,0.75) quantile(tmp,0.025) quantile(tmp,0.975)};
            for j = 1:6                
                FSW_Cases{ind,j} = tmp{j};
            end
            
            tmp = NewCases_M./(TotM)*1000;
            tmp = {year mean(tmp) quantile(tmp,0.25) quantile(tmp,0.75) quantile(tmp,0.025) quantile(tmp,0.975)};
            for j = 1:6                
                M_Cases{ind,j} = tmp{j};
            end
            
            NewCases_FSW = 0;
            NewCases_M = 0;
        end 
       
    end

    
    PlotForecast(Res,record,go_back)
    SavePath = 'H:\My Documents\PhD Work\Imperial Project\Meetings\Images\';

    saveas(gcf,[SavePath Names{i} '_ForecPeter50years.eps'], 'psc2')

    disp(i)
    disp(Names{i})
    
    SavePath = 'H:\My Documents\PhD Work\Imperial Project\Meetings\xlsFiles';
    xlswrite([SavePath Names{i} '_PrevalenceEstFSW.xls'], FSW_Prev, 'Prevalence Estimates FSW', 'A1');
    xlswrite([SavePath Names{i} '_PrevalenceEstClients.xls'], M_Prev, 'Prevalence Estimates M', 'A1');
    xlswrite([SavePath Names{i} '_NbYearlyCasesPer1000_FSW.xls'], FSW_Cases, 'New cases per 1000 FSW', 'A1');
    xlswrite([SavePath Names{i} '_NbYearlyCasesPer1000_Clients.xls'], M_Prev, 'New cases per 1000 clients', 'A1');

end








