function [] = PlotParsHIV(Res,FtPath)

Data = Res.Data;
Parameters = Res.Parameters;
Names = Parameters.Names.Estimated;
k = length(Parameters.Names.All)-4;

figure(1)
indplot = 0;
for i = 1:length(Parameters.Names.All)
    if not(strcmp(Names{i}(end-1:end),'Ft'))
        indplot = indplot+1;
        subplot(ceil(sqrt(k)),ceil(sqrt(k)),indplot)
        if Parameters.(Names{i}).Estimated;
            ind = Parameters.(Names{i}).Index;
            minval = min(Res.Thetas(ind,:));
            maxval = max(Res.Thetas(ind,:));
            TheoMin = Parameters.(Names{i}).Min;
            TheoMax = Parameters.(Names{i}).Max;
            TheoX = randn(10000,1)*(TheoMax-TheoMin)/4 + (TheoMax+TheoMin)/2;
            xis = minval:(maxval-minval)/1000:maxval;
            [Theofi,xis]=ksdensity(TheoX,xis); 
            [Empfi,xis]=ksdensity(Res.Thetas(ind,:),xis); 
            plot(xis,Theofi,'g')
            hold on
            plot(xis,Empfi)
            hold off
            title(Names{i})
        else
            ind = Parameters.(Names{i}).Index;
            minval = Parameters.(Names{i}).Min;
            maxval = Parameters.(Names{i}).Max;
            TheoMin = Parameters.(Names{i}).Min;
            TheoMax = Parameters.(Names{i}).Max;
            TheoX = randn(10000,1)*(TheoMax-TheoMin)/4 + (TheoMax+TheoMin)/2;
            xis = minval:(maxval-minval)/1000:maxval;
            [Theofi,xis]=ksdensity(TheoX,xis); 
            plot(xis,Theofi,'g')
            hold on
            plot(Parameters.(Names{i}).Value*ones(2,1),[0 max(Theofi)])
            hold off
            title(Names{i})
        end
    end
end
figure(2)
indplot = 0;
for i = 1:length(Parameters.Names.All)
    if not(strcmp(Names{i}(end-1:end),'Ft'))
        indplot = indplot+1;
        subplot(ceil(sqrt(k)),ceil(sqrt(k)),indplot)
        if Parameters.(Names{i}).Estimated;
            ind = Parameters.(Names{i}).Index;
            minval = min(Res.Thetas(ind,:));
            maxval = max(Res.Thetas(ind,:));
            TheoMin = Parameters.(Names{i}).Min;
            TheoMax = Parameters.(Names{i}).Max;
            TheoX = randn(10000,1)*(TheoMax-TheoMin)/4 + (TheoMax+TheoMin)/2;
            xis = minval:(maxval-minval)/1000:maxval;
            plot(Res.Thetas(ind,:))
            title(Names{i})
        else
            ind = Parameters.(Names{i}).Index;
            minval = Parameters.(Names{i}).Min;
            maxval = Parameters.(Names{i}).Max;
            TheoMin = Parameters.(Names{i}).Min;
            TheoMax = Parameters.(Names{i}).Max;
            TheoX = randn(10000,1)*(TheoMax-TheoMin)/4 + (TheoMax+TheoMin)/2;
            xis = minval:(maxval-minval)/1000:maxval;
            [Theofi,xis]=ksdensity(TheoX,xis); 
             plot(xis,Theofi,'g')
            hold on
            plot(Parameters.(Names{i}).Value*ones(2,1),[0 max(Theofi)])
            hold off
            title(Names{i})
        end
    end
end