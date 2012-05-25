function ind = PickRandComp( x , DensityModel, Parameters)

NComps = DensityModel.NComponents;

liks = [];
for i = 1:NComps
%     liks(i) = DensityModel.PComponents(i);
    liks(i) = DensityModel.PComponents(i)*(mvnpdf(x,DensityModel.mu(i,:)',squeeze(DensityModel.Sigma(:,:,i))));
%     liks(i) = (mvnpdf(x,DensityModel.mu(i,:)',squeeze(DensityModel.Sigma(:,:,i))));
end

NoPlot = 1;
if length(x) == 2
    xmin = [0.3 0.3*10^-5];
    xmax = [75  1.8*10^-5];
    NoPlot = 0;
elseif length(x) == 3
    xmin = [0.3 0.8 0.3*10^-5];
    xmax = [75  1.2 1.8*10^-5];
    NoPlot = 0;
end

if not(NoPlot)
    ScaledLogXmin = Scale(log(xmin'),Parameters);
    ScaledLogXmax = Scale(log(xmax'),Parameters);

    for i = 1:length(x)
        subplot(6*length(x),1,(i-1)*6+1:i*6-1)
        var = i;
        plot(x(var),0,'og')
        hold on
        xs = ScaledLogXmin(var):0.01:ScaledLogXmax(var);
        for j = 1:NComps
            temp = DensityModel.PComponents(j)*normpdf(xs,DensityModel.mu(j,var),sqrt(DensityModel.Sigma(var,var,j)));
            plot(xs,temp)
        end
        hold off
        subplot(6*length(x),1,i*6)
        plot(exp(x(var)),0,'og')
        xlim([xmin(i) xmax(i)])
        pause(0.01)
    end
end

ind = PickRandInd(liks);