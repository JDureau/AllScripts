function ResStats = PlotPostForqNew(Mode,q,s,Delta)


files = {'PostThreeMethods.mat','PostThreeMethodsAffine.mat','PostThreeMethodsStep.mat'};
SavePath = '/Users/dureaujoseph/Documents/Taf/These/Matlab Scripts/AllData/Avahan/';


Stats = {'ampl aft 2003','asympt'};

NbBS = 100;

ResTmp = {};

ResTmp{1,2} = 'dBRSim';
ResTmp{1,3} = 'dBRSim';
ResTmp{1,4} = 'dBRSim';
ResTmp{1,5} = 'dBRSim';
ResTmp{1,6} = 'dBRSim';
% ResTmp{1,7} = 'dBRSim';
% ResTmp{1,7} = 'AffSim';
% ResTmp{1,8} = 'AffSim';
% ResTmp{1,9} = 'AffSim';
% ResTmp{1,10} = 'StepSim';
% ResTmp{1,11} = 'StepSim';
% ResTmp{1,12} = 'StepSim';
% ResTmp{1,13} = 'StepSim';
ResTmp{2,2} = 'dBRMod';
ResTmp{2,3} = 'sBRMod';
ResTmp{2,4} = 'BMMod';
ResTmp{2,5} = 'SigmMod';
ResTmp{2,6} = 'sSigmMod';
% ResTmp{2,6} = 'dBRMod';
% ResTmp{2,7} = 'sBRMod';
% ResTmp{2,8} = 'BMMod';
% ResTmp{2,9} = 'SigmMod';
% ResTmp{2,10} = 'dBRMod';
% ResTmp{2,11} = 'sBRMod';
% ResTmp{2,12} = 'BMMod';
% ResTmp{2,13} = 'SigmMod';
ResTmp{3,1} = 'Bias' ;
ResTmp{4,1} = 'Variance' ;
ResTmp{5,1} = 'MSE' ;
ResTmp{6,1} = 'AUC' ;
ResTmp{7,1} = 'TPR' ;
ResTmp{8,1} = 'FPR' ;


ResDelta2003 = ResTmp;
ResCU2010 = ResTmp;
ResDelta2003{1,1} = ['Delta=' num2str(Delta)];
ResCU2010{1,1} = ['CU2010=' num2str(Delta)];
Pts = {};
ROCs = [];


for f = 1:1%3
    f
    load([SavePath files{f}])
    
    % Stats = {'ampl','ampl. aft. 1995','ampl aft 2003','asympt'};
    Tabs = {};


    inds = 1:size(Res.ampls,2);


    if strcmp(Mode,'ROC')
        indq = find(Res.qs == q);
        Tabs{1} = Res.amplsROC(:,indq,inds);
        Tabs{2} = Res.partampls1995ROC(:,indq,inds);
        Tabs{3} = Res.partampls2003ROC(:,indq,inds);
        Tabs{4} = Res.asymptsROC(:,indq,inds);
    elseif strcmp(Mode,'Post')
        Tabs{1} = Res.amplsPost(:,inds);
        Tabs{2} = Res.partampls1995Post(:,inds);
        Tabs{3} = Res.partampls2003Post(:,inds);
        Tabs{4} = Res.asymptsPost(:,inds);
    else
        Tabs{1} = Res.ampls(:,inds);
        Tabs{2} = Res.partampls1995(:,inds);
        Tabs{3} = Res.partampls2003(:,inds);
        Tabs{4} = Res.asympts(:,inds);
    end




    ParsBiases = [];
    ParsVars = [];
    ParsMSEs = [];

    Corrs = [];
    Biases = [];
    MSEs = [];
    Vars = [];
    inds = [1 3 2 0 4 5];
    for i = 3:4
        for j = [1 2 3 5 6]
            Pts{f,i,inds(j)} = [];
            for kbs = 1:NbBS
                bsinds = ceil(rand(1,size(Tabs{i},3))*size(Tabs{i},3));
                
                Corrs(inds(j),i,kbs) = corr(Tabs{i}(4,bsinds)',Tabs{i}(j,bsinds)');
                Biases(inds(j),i,kbs) = mean(Tabs{i}(j,bsinds)-Tabs{i}(4,bsinds)); 
                MSEs(inds(j),i,kbs) = mean((Tabs{i}(4,bsinds)-Tabs{i}(j,bsinds)).^2);
                Vars(inds(j),i,kbs) = MSEs(inds(j),i,kbs)-Biases(inds(j),i,kbs).^2;

                if i == 3
                    quan = Delta;
                elseif i == 4
                    quan = Delta;
                else
                    quan = 0;
                end
                indsinf = find(Tabs{i}(4,bsinds)<=quan);
                indssup = find(Tabs{i}(4,bsinds)>quan);
                xs = [];
                ys = [];
                all = [Tabs{i}(1,:) Tabs{i}(2,:) Tabs{i}(3,:)];
                m = min(all);
                M = max(all);
                deltap = sort(all);
                deltap = (deltap(1:end-1)+deltap(2:end))/2;
                deltap = [(min(all)-mean(diff(deltap))) deltap (max(all)+mean(diff(deltap)))];
                for k = 1:length(deltap)
                    xs(k) = sum(Tabs{i}(j,bsinds(indsinf))>deltap(k))/length(indsinf);
                    ys(k) = sum(Tabs{i}(j,bsinds(indssup))>deltap(k))/length(indssup);
                end
                [b,indk] = min(abs(deltap-quan));
                TruePosRate(inds(j),i,kbs) = ys(indk);
                FalsePosRate(inds(j),i,kbs) = xs(indk);
                
                
                ROCs(f,inds(j),i,kbs) = 0;
                for k = 1:length(deltap)-1
                    ROCs(f,inds(j),i,kbs) = ROCs(f,inds(j),i,kbs) + (xs(k)-xs(k+1))*(ys(k+1)+ys(k))/2;
                end
%                 ROCs(inds(j),i,kbs) = ROCs(inds(j),i,kbs)-0.5;
                Pts{f,i,inds(j)}(kbs,:,:)=[xs;ys];

            end
        end
    end
    
%     i=3
%     j1=1
%     j2=2
%     plot(squeeze(Tabs{i}(j1,1,:)-Tabs{i}(4,1,:)),squeeze(Tabs{i}(j2,1,:)-Tabs{i}(4,1,:)),'.')

    % Find quantile that enables 5% false positives and give corresponding true
    % positives rate
%     Rightqs = [];
%     qs = Res.qs;
%     FalsePosRate = [];
%     TruePosRate = [];
%     TempFalsePosRate = [];
%     TempTruePosRate = [];
%     qMSEs = [];
%     qBiases = [];
%     nbyes = 0;
%     nbno = 0;
%     for i = 3:4
%         for j = 1:3
%             TruePosRate(inds(j),i)  = 0;
%             FalsePosRate(inds(j),i)  = 0;
%             nbyes = 0;
%             nbno = 0;
%             for kbs = 1:NbBS
% 
%             
%                 
%                
%                 if i == 3
%                     quan = Delta;
%                 elseif i == 4
%                     quan = Asympt;
%                 else
%                     quan = 0;
%                 end
% 
% %                 test = 0;
% %                 while not(test)
% %                     if Tabs{i}(4,bsind)>quan
% %                         test = 1;
% %                     end
% %                 end
%                 
%                 bsind = ceil(rand(1,1)*size(Tabs{i},3));
%                 bsothers = [1:bsind-1 bsind+1:size(Tabs{i},3)];
%     
%                 indsinf = find(Tabs{i}(4,bsothers)<=quan);
%                 indssup = find(Tabs{i}(4,bsothers)>quan);
% 
% 
%                 for k = 1:length(qs)
%                     inds = 1:size(Res.ampls,2);
%                     Tabs{1} = Res.amplsROC(:,k,inds);
%                     Tabs{2} = Res.partampls1995ROC(:,k,inds);
%                     Tabs{3} = Res.partampls2003ROC(:,k,inds);
%                     Tabs{4} = Res.asymptsROC(:,k,inds);
% 
%                     xs = [];
%                     ys = [];
%                     all = [Tabs{i}(1,:) Tabs{i}(2,:) Tabs{i}(3,:)];
%                     m = min(all);
%                     M = max(all);
%             %         deltap = sort((all(1:end-1)+all(2:end))/2);
%                     deltap = sort(all);
%                     deltap = (deltap(1:end-1)+deltap(2:end))/2;
%                     deltap = [(min(all)-mean(diff(deltap))) deltap (max(all)+mean(diff(deltap)))];
% 
%                     TempFalsePosRate(k) = sum(Tabs{i}(j,bsothers(indsinf))>quan)/length(bsothers(indsinf));
%                     TempTruePosRate(k) = sum(Tabs{i}(j,bsothers(indssup))>quan)/length(bsothers(indssup));
%                 end
%                 [b,ind] = min(abs(TempFalsePosRate-0.05)+5*(TempFalsePosRate==0));
%         %         ind = 5;
%                 Tabs{1} = Res.amplsROC(:,ind,inds);
%                 Tabs{2} = Res.partampls1995ROC(:,ind,inds);
%                 Tabs{3} = Res.partampls2003ROC(:,ind,inds);
%                 Tabs{4} = Res.asymptsROC(:,ind,inds);
%                 inds = [1 3 2];
%                 if Tabs{i}(4,bsind)>quan
%                     nbyes = nbyes + 1;
%                     if (Tabs{i}(j,bsind)>quan)&&(Tabs{i}(4,bsind)>quan)
%                         TruePosRate(inds(j),i) = TruePosRate(inds(j),i) + 1;
%                     end
%                 else
%                     nbno = nbno + 1;
%                     if (Tabs{i}(j,bsind)>quan)&&(Tabs{i}(4,bsind)<=quan)
%                         FalsePosRate(inds(j),i) = FalsePosRate(inds(j),i) + 1;
%                     end
%                 end
%                 Rightqs(inds(j),i,kbs) = qs(ind);
% %                 qMSEs(inds(j),i) = mean((Tabs{i}(4,:)-Tabs{i}(j,:)).^2);
% %                 qBiases(inds(j),i) = mean((Tabs{i}(j,:)-Tabs{i}(4,:)));
%             end
%             TruePosRate(inds(j),i)  = TruePosRate(inds(j),i)/nbyes;
%             FalsePosRate(inds(j),i)  = FalsePosRate(inds(j),i)/nbno;
%         end
%     end


    Pars = Res.ParsMeds;
    Pars2p5 = Res.Pars2p5;
    Pars97p5 = Res.Pars97p5;
%     Pars = Res.Pars;
    Names = Res.Parameters.Names.Estimated;
    propin = [];
    for i = 1:size(Pars,2)
        for j = 1:3
            inddet = Res.Parameters.(Names{i}).Index;
            
            ParsBiases(inds(j),inddet) = mean(Pars(4,inddet,:)-Pars(j,inddet,:));
            ParsVars(inds(j),inddet) = mean((Pars(4,inddet,:)-Pars(j,inddet,:)).^2);
            ParsMSEs(inds(j),inddet) = ParsVars(inds(j),inddet) + (ParsBiases(inds(j),inddet))^2;
            propin(inds(j),inddet) = mean((squeeze((Pars(4,inddet,:)))>squeeze(Pars2p5(j,inddet,:))).*(squeeze(Pars(4,inddet,:))<squeeze(Pars97p5(j,inddet,:))));
        end
    end
%     for i = 1:size(Pars,2)
%         inddet = Res.Parameters.(Names{i}).Index;
%         disp(Names{i});
%         ParsBiases(:,inddet)
%     %     ParsVars(:,inddet)
%         ParsMSEs(:,inddet)
%     end
% propin(:,13)
% ParsBiases(:,Res.Parameters.BRmm1.Index)


    % subplot(4,2,1:2)
    % bar(Biases(:,3:4)')
    % set(gca,'XTickLabel',Stats)
    % ylim([-0.4 0.2])
    % ylabel('Biases')
    % Biases(:,[1 3 4])
    % 
    % subplot(4,2,3:4)
    % bar(Vars(:,3:4)')
    % set(gca,'XTickLabel',Stats)
    % ylim([0 0.1])
    % ylabel('Vars')
    % Biases(:,[1 3 4])
    % 
    % subplot(4,2,5:6)
    % bar((MSEs(:,3:4))')
    % set(gca,'XTickLabel',Stats)
    % ylim([0 0.1])
    % ylabel('MSEs')
    % legend('Det. BR','Sto. BR',' BM')
    % MSEs(:,[1 3 4])
    % 
    % 
    % subplot(4,2,7:8)
    % bar(ROCs(:,3:4)')
    % ylim([0 0.5])
    % set(gca,'XTickLabel',Stats)
    % ylabel('ROCs-0.5')


%     subplot(4,2,1:2)
%     bar(sqrt(MSEs(:,3:4))')
%     set(gca,'XTickLabel',Stats)
%     ylim([0 0.4])
%     ylabel('sqrt(MSEs)')
% 
% 
%     subplot(4,2,3:4)
%     bar(ROCs(:,3:4)')
%     ylim([0 0.5])
%     set(gca,'XTickLabel',Stats)
%     ylabel('ROCs-0.5')
% 
%     subplot(4,2,5:6)
%     bar(sqrt(qMSEs(:,3:4))')
%     ylim([0 0.4])
%     set(gca,'XTickLabel',Stats)
%     ylabel('sqrt(qMSEs)')
%     % bar(FalsePosRate(:,3:4)')
%     % ylim([0 1])
%     % set(gca,'XTickLabel',Stats)
%     % ylabel('False pos. rate')
%     title([num2str(Rightqs(:,3)') '                  ' num2str(Rightqs(:,4)')])
% 
% 
%     subplot(4,2,7:8)
%     bar(TruePosRate(:,3:4)')
%     ylim([0 1])
%     set(gca,'XTickLabel',Stats)
%     ylabel('True pos. rate')
%     title([num2str(Rightqs(:,3)') '                  ' num2str(Rightqs(:,4)')])

%     i = 4;
%     % disp('Bias:')
%     % Biases(:,i)'
%     % disp('Vars:')
%     % Vars(:,i)'
%     disp('MSEs:')
%     MSEs(:,i)'
%     disp('ROCs:')
%     ROCs(:,i)'
%     disp('qs:')
%     Rightqs(:,i)'
%     % disp('FalsePos:')
%     % FalsePosRate(:,i)'
%     % disp('qMSEs:')
%     % qMSEs(:,i)'
%     disp('TruePos:')
%     TruePosRate(:,i)'


    ResDelta2003{3,2+4*(f-1)}=[num2str(mean(Biases(1,3,:)),2) ' [' num2str(quantile(Biases(1,3,:),0.025),2) ';' num2str(quantile(Biases(1,3,:),0.975),2) ']'];
    ResDelta2003{3,3+4*(f-1)}=[num2str(mean(Biases(2,3,:)),2) ' [' num2str(quantile(Biases(2,3,:),0.025),2) ';' num2str(quantile(Biases(2,3,:),0.975),2) ']'];
    ResDelta2003{3,4+4*(f-1)}=[num2str(mean(Biases(3,3,:)),2) ' [' num2str(quantile(Biases(3,3,:),0.025),2) ';' num2str(quantile(Biases(3,3,:),0.975),2) ']'];
    ResDelta2003{3,5+4*(f-1)}=[num2str(mean(Biases(4,3,:)),2) ' [' num2str(quantile(Biases(4,3,:),0.025),2) ';' num2str(quantile(Biases(4,3,:),0.975),2) ']'];
    ResDelta2003{3,6+4*(f-1)}=[num2str(mean(Biases(5,3,:)),2) ' [' num2str(quantile(Biases(5,3,:),0.025),2) ';' num2str(quantile(Biases(5,3,:),0.975),2) ']'];
    ResDelta2003{4,2+4*(f-1)}=[num2str(mean(Vars(1,3,:)),2) ' [' num2str(quantile(Vars(1,3,:),0.025),2) ';' num2str(quantile(Vars(1,3,:),0.975),2) ']'];
    ResDelta2003{4,3+4*(f-1)}=[num2str(mean(Vars(2,3,:)),2) ' [' num2str(quantile(Vars(2,3,:),0.025),2) ';' num2str(quantile(Vars(2,3,:),0.975),2) ']'];
    ResDelta2003{4,4+4*(f-1)}=[num2str(mean(Vars(3,3,:)),2) ' [' num2str(quantile(Vars(3,3,:),0.025),2) ';' num2str(quantile(Vars(3,3,:),0.975),2) ']'];
    ResDelta2003{4,5+4*(f-1)}=[num2str(mean(Vars(4,3,:)),2) ' [' num2str(quantile(Vars(4,3,:),0.025),2) ';' num2str(quantile(Vars(4,3,:),0.975),2) ']'];
    ResDelta2003{4,6+4*(f-1)}=[num2str(mean(Vars(5,3,:)),2) ' [' num2str(quantile(Vars(5,3,:),0.025),2) ';' num2str(quantile(Vars(5,3,:),0.975),2) ']'];
    ResDelta2003{5,2+4*(f-1)}=[num2str(mean(MSEs(1,3,:)),2) ' [' num2str(quantile(MSEs(1,3,:),0.025),2) ';' num2str(quantile(MSEs(1,3,:),0.975),2) ']'];
    ResDelta2003{5,3+4*(f-1)}=[num2str(mean(MSEs(2,3,:)),2) ' [' num2str(quantile(MSEs(2,3,:),0.025),2) ';' num2str(quantile(MSEs(2,3,:),0.975),2) ']'];
    ResDelta2003{5,4+4*(f-1)}=[num2str(mean(MSEs(3,3,:)),2) ' [' num2str(quantile(MSEs(3,3,:),0.025),2) ';' num2str(quantile(MSEs(3,3,:),0.975),2) ']'];
    ResDelta2003{5,5+4*(f-1)}=[num2str(mean(MSEs(4,3,:)),2) ' [' num2str(quantile(MSEs(4,3,:),0.025),2) ';' num2str(quantile(MSEs(4,3,:),0.975),2) ']'];
    ResDelta2003{5,6+4*(f-1)}=[num2str(mean(MSEs(5,3,:)),2) ' [' num2str(quantile(MSEs(5,3,:),0.025),2) ';' num2str(quantile(MSEs(5,3,:),0.975),2) ']'];
    ResDelta2003{6,2+4*(f-1)}=[num2str(mean(ROCs(f,1,3,:)),2) ' [' num2str(quantile(ROCs(f,1,3,:),0.025),2) ';' num2str(quantile(ROCs(f,1,3,:),0.975),2) ']'];
    ResDelta2003{6,3+4*(f-1)}=[num2str(mean(ROCs(f,2,3,:)),2) ' [' num2str(quantile(ROCs(f,2,3,:),0.025),2) ';' num2str(quantile(ROCs(f,2,3,:),0.975),2) ']'];
    ResDelta2003{6,4+4*(f-1)}=[num2str(mean(ROCs(f,3,3,:)),2) ' [' num2str(quantile(ROCs(f,3,3,:),0.025),2) ';' num2str(quantile(ROCs(f,3,3,:),0.975),2) ']'];
    ResDelta2003{6,5+4*(f-1)}=[num2str(mean(ROCs(f,4,3,:)),2) ' [' num2str(quantile(ROCs(f,4,3,:),0.025),2) ';' num2str(quantile(ROCs(f,4,3,:),0.975),2) ']'];
    ResDelta2003{6,6+4*(f-1)}=[num2str(mean(ROCs(f,4,3,:)),2) ' [' num2str(quantile(ROCs(f,5,3,:),0.025),2) ';' num2str(quantile(ROCs(f,5,3,:),0.975),2) ']'];
    ResDelta2003{7,2+4*(f-1)}=[num2str(mean(TruePosRate(1,3,:)),2) ' [' num2str(quantile(TruePosRate(1,3,:),0.025),2) ';' num2str(quantile(TruePosRate(1,3,:),0.975),2) ']'];
    ResDelta2003{7,3+4*(f-1)}=[num2str(mean(TruePosRate(2,3,:)),2) ' [' num2str(quantile(TruePosRate(2,3,:),0.025),2) ';' num2str(quantile(TruePosRate(2,3,:),0.975),2) ']'];
    ResDelta2003{7,4+4*(f-1)}=[num2str(mean(TruePosRate(3,3,:)),2) ' [' num2str(quantile(TruePosRate(3,3,:),0.025),2) ';' num2str(quantile(TruePosRate(3,3,:),0.975),2) ']'];
    ResDelta2003{7,5+4*(f-1)}=[num2str(mean(TruePosRate(4,3,:)),2) ' [' num2str(quantile(TruePosRate(4,3,:),0.025),2) ';' num2str(quantile(TruePosRate(4,3,:),0.975),2) ']'];
    ResDelta2003{7,6+4*(f-1)}=[num2str(mean(TruePosRate(5,3,:)),2) ' [' num2str(quantile(TruePosRate(5,3,:),0.025),2) ';' num2str(quantile(TruePosRate(5,3,:),0.975),2) ']'];
    ResDelta2003{8,2+4*(f-1)}=[num2str(mean(FalsePosRate(1,3,:)),2) ' [' num2str(quantile(FalsePosRate(1,3,:),0.025),2) ';' num2str(quantile(FalsePosRate(1,3,:),0.975),2) ']'];
    ResDelta2003{8,3+4*(f-1)}=[num2str(mean(FalsePosRate(2,3,:)),2) ' [' num2str(quantile(FalsePosRate(2,3,:),0.025),2) ';' num2str(quantile(FalsePosRate(2,3,:),0.975),2) ']'];
    ResDelta2003{8,4+4*(f-1)}=[num2str(mean(FalsePosRate(3,3,:)),2) ' [' num2str(quantile(FalsePosRate(3,3,:),0.025),2) ';' num2str(quantile(FalsePosRate(3,3,:),0.975),2) ']'];
    ResDelta2003{8,5+4*(f-1)}=[num2str(mean(FalsePosRate(4,3,:)),2) ' [' num2str(quantile(FalsePosRate(4,3,:),0.025),2) ';' num2str(quantile(FalsePosRate(4,3,:),0.975),2) ']'];
    ResDelta2003{8,6+4*(f-1)}=[num2str(mean(FalsePosRate(5,3,:)),2) ' [' num2str(quantile(FalsePosRate(5,3,:),0.025),2) ';' num2str(quantile(FalsePosRate(5,3,:),0.975),2) ']'];
    
    ResCU2010{3,2+4*(f-1)}=[num2str(mean(Biases(1,4,:)),2) ' [' num2str(quantile(Biases(1,4,:),0.025),2) ';' num2str(quantile(Biases(1,4,:),0.975),2) ']'];
    ResCU2010{3,3+4*(f-1)}=[num2str(mean(Biases(2,4,:)),2) ' [' num2str(quantile(Biases(2,4,:),0.025),2) ';' num2str(quantile(Biases(2,4,:),0.975),2) ']'];
    ResCU2010{3,4+4*(f-1)}=[num2str(mean(Biases(3,4,:)),2) ' [' num2str(quantile(Biases(3,4,:),0.025),2) ';' num2str(quantile(Biases(3,4,:),0.975),2) ']'];
    ResCU2010{3,5+4*(f-1)}=[num2str(mean(Biases(4,4,:)),2) ' [' num2str(quantile(Biases(4,4,:),0.025),2) ';' num2str(quantile(Biases(4,4,:),0.975),2) ']'];
    ResCU2010{3,6+4*(f-1)}=[num2str(mean(Biases(5,4,:)),2) ' [' num2str(quantile(Biases(5,4,:),0.025),2) ';' num2str(quantile(Biases(5,4,:),0.975),2) ']'];
    ResCU2010{4,2+4*(f-1)}=[num2str(mean(Vars(1,4,:)),2) ' [' num2str(quantile(Vars(1,4,:),0.025),2) ';' num2str(quantile(Vars(1,4,:),0.975),2) ']'];
    ResCU2010{4,3+4*(f-1)}=[num2str(mean(Vars(2,4,:)),2) ' [' num2str(quantile(Vars(2,4,:),0.025),2) ';' num2str(quantile(Vars(2,4,:),0.975),2) ']'];
    ResCU2010{4,4+4*(f-1)}=[num2str(mean(Vars(3,4,:)),2) ' [' num2str(quantile(Vars(3,4,:),0.025),2) ';' num2str(quantile(Vars(3,4,:),0.975),2) ']'];
    ResCU2010{4,5+4*(f-1)}=[num2str(mean(Vars(4,4,:)),2) ' [' num2str(quantile(Vars(4,4,:),0.025),2) ';' num2str(quantile(Vars(4,4,:),0.975),2) ']'];
    ResCU2010{4,6+4*(f-1)}=[num2str(mean(Vars(5,4,:)),2) ' [' num2str(quantile(Vars(5,4,:),0.025),2) ';' num2str(quantile(Vars(5,4,:),0.975),2) ']'];
    ResCU2010{5,2+4*(f-1)}=[num2str(mean(MSEs(1,4,:)),2) ' [' num2str(quantile(MSEs(1,4,:),0.025),2) ';' num2str(quantile(MSEs(1,4,:),0.975),2) ']'];
    ResCU2010{5,3+4*(f-1)}=[num2str(mean(MSEs(2,4,:)),2) ' [' num2str(quantile(MSEs(2,4,:),0.025),2) ';' num2str(quantile(MSEs(2,4,:),0.975),2) ']'];
    ResCU2010{5,4+4*(f-1)}=[num2str(mean(MSEs(3,4,:)),2) ' [' num2str(quantile(MSEs(3,4,:),0.025),2) ';' num2str(quantile(MSEs(3,4,:),0.975),2) ']'];
    ResCU2010{5,5+4*(f-1)}=[num2str(mean(MSEs(4,4,:)),2) ' [' num2str(quantile(MSEs(4,4,:),0.025),2) ';' num2str(quantile(MSEs(4,4,:),0.975),2) ']'];
    ResCU2010{5,6+4*(f-1)}=[num2str(mean(MSEs(5,4,:)),2) ' [' num2str(quantile(MSEs(5,4,:),0.025),2) ';' num2str(quantile(MSEs(5,4,:),0.975),2) ']'];
    ResCU2010{6,2+4*(f-1)}=[num2str(mean(ROCs(f,1,4,:)),2) ' [' num2str(quantile(ROCs(f,1,4,:),0.025),2) ';' num2str(quantile(ROCs(f,1,4,:),0.975),2) ']'];
    ResCU2010{6,3+4*(f-1)}=[num2str(mean(ROCs(f,2,4,:)),2) ' [' num2str(quantile(ROCs(f,2,4,:),0.025),2) ';' num2str(quantile(ROCs(f,2,4,:),0.975),2) ']'];
    ResCU2010{6,4+4*(f-1)}=[num2str(mean(ROCs(f,3,4,:)),2) ' [' num2str(quantile(ROCs(f,3,4,:),0.025),2) ';' num2str(quantile(ROCs(f,3,4,:),0.975),2) ']'];
    ResCU2010{6,5+4*(f-1)}=[num2str(mean(ROCs(f,4,4,:)),2) ' [' num2str(quantile(ROCs(f,4,4,:),0.025),2) ';' num2str(quantile(ROCs(f,4,4,:),0.975),2) ']'];
    ResCU2010{7,2+4*(f-1)}=[num2str(mean(TruePosRate(1,4,:)),2) ' [' num2str(quantile(TruePosRate(1,4,:),0.025),2) ';' num2str(quantile(TruePosRate(1,4,:),0.975),2) ']'];
    ResCU2010{7,3+4*(f-1)}=[num2str(mean(TruePosRate(2,4,:)),2) ' [' num2str(quantile(TruePosRate(2,4,:),0.025),2) ';' num2str(quantile(TruePosRate(2,4,:),0.975),2) ']'];
    ResCU2010{7,4+4*(f-1)}=[num2str(mean(TruePosRate(3,4,:)),2) ' [' num2str(quantile(TruePosRate(3,4,:),0.025),2) ';' num2str(quantile(TruePosRate(3,4,:),0.975),2) ']'];
    ResCU2010{7,5+4*(f-1)}=[num2str(mean(TruePosRate(4,4,:)),2) ' [' num2str(quantile(TruePosRate(4,4,:),0.025),2) ';' num2str(quantile(TruePosRate(4,4,:),0.975),2) ']'];
    ResCU2010{7,6+4*(f-1)}=[num2str(mean(TruePosRate(5,4,:)),2) ' [' num2str(quantile(TruePosRate(5,4,:),0.025),2) ';' num2str(quantile(TruePosRate(5,4,:),0.975),2) ']'];
    ResCU2010{8,2+4*(f-1)}=[num2str(mean(FalsePosRate(1,4,:)),2) ' [' num2str(quantile(FalsePosRate(1,4,:),0.025),2) ';' num2str(quantile(FalsePosRate(1,4,:),0.975),2) ']'];
    ResCU2010{8,3+4*(f-1)}=[num2str(mean(FalsePosRate(2,4,:)),2) ' [' num2str(quantile(FalsePosRate(2,4,:),0.025),2) ';' num2str(quantile(FalsePosRate(2,4,:),0.975),2) ']'];
    ResCU2010{8,4+4*(f-1)}=[num2str(mean(FalsePosRate(3,4,:)),2) ' [' num2str(quantile(FalsePosRate(3,4,:),0.025),2) ';' num2str(quantile(FalsePosRate(3,4,:),0.975),2) ']'];
    ResCU2010{8,5+4*(f-1)}=[num2str(mean(FalsePosRate(4,4,:)),2) ' [' num2str(quantile(FalsePosRate(4,4,:),0.025),2) ';' num2str(quantile(FalsePosRate(4,4,:),0.975),2) ']'];
    ResCU2010{8,6+4*(f-1)}=[num2str(mean(FalsePosRate(5,4,:)),2) ' [' num2str(quantile(FalsePosRate(5,4,:),0.025),2) ';' num2str(quantile(FalsePosRate(5,4,:),0.975),2) ']'];
    
end

AvahanFile = '/Users/dureaujoseph/Dropbox/Taf/Notes/SecondPaper/';

if strcmp(s,'Delta2003')
    cell2csv([AvahanFile 'RunsOutputsDelta2003_q_' num2str(q) '__Delta' num2str(Delta) '.csv'], ResDelta2003)
elseif strcmp(s,'CU2010')
    cell2csv([AvahanFile 'RunsOutputsCU2010_q_' num2str(q) '__Delta' num2str(Delta) '.csv'], ResCU2010)    
end


    

% disp('qBiases:')
% qBiases(:,i)'
% 
% subplot(6,4,1:4)
% bar(Corrs')
% set(gca,'XTickLabel',Stats)
% ylim([0 1])
% ylabel('Correlation')
% 
% subplot(6,4,5:8)
% bar(Biases')
% set(gca,'XTickLabel',Stats)
% ylim([-0.4 0.2])
% ylabel('Biases')
% Biases(:,[1 3 4])
% 
% subplot(6,4,9:12)
% bar((MSEs)')
% set(gca,'XTickLabel',Stats)
% ylim([0 0.2])
% ylabel('MSEs')
% legend('Det. BR','Sto. BR',' BM')
% MSEs(:,[1 3 4])
% 
% subplot(6,4,13:16)
% bar(Vars')
% set(gca,'XTickLabel',Stats)
% ylabel('Vars')
% 
% subplot(6,4,21:24)
% bar(ROCs')
% ylim([0 0.5])
% set(gca,'XTickLabel',Stats)
% ylabel('ROCs-0.5')
% ROCs(:,[1 3 4])
% 
% for i = 1:4
%    subplot(6,4,16+i) 
%    plot(Pts{i,1}(1,:),Pts{i,1}(2,:),'b')
%    hold on
%    plot(Pts{i,2}(1,:),Pts{i,2}(2,:),'g')
%    plot(Pts{i,3}(1,:),Pts{i,3}(2,:),'r')
%    hold off
% end

ResStats = struct();
ResStats.Corrs = Corrs;
ResStats.Biases = Biases;
ResStats.MSEs = MSEs;
ResStats.Vars = Vars;
ResStats.ROCs = ROCs;
ResStats.quant = quan;
ResStats.ParsBiases = ParsBiases;
ResStats.ParsVars = ParsVars;
ResStats.ParsMSEs = ParsMSEs;
ResStats.Pts = Pts;