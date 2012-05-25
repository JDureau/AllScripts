function [] = PlotTracePlots(Res,NbSamples,MaxLag)

% for autocorrelation
if nargin == 1
    NbSamples = size(Res.Thetas,2);
    MaxLag = min(500,NbSamples);
end

dim = Res.Parameters.NbParsEstimated;
figure(1)
clf
if dim <6
    for i = 1:dim
        for j = i+1:dim
            subplot(dim,dim,(i-1)*dim+j)
            plot(Res.Thetas(j,:),Res.Thetas(i,:),'.')
            if i == 1
                title(Res.Parameters.Names.Estimated{j})
            end
            subplot(dim,dim,(j-1)*dim+i)
            plot(Res.TransfThetas(i,:),Res.TransfThetas(j,:),'.')
            if i == 1
                ylabel(Res.Parameters.Names.Estimated{j})
            end
        end
        subplot(dim,dim,(i-1)*dim+i)
        plot(Res.Thetas(i,:))
        if i == 1
            ylabel(Res.Parameters.Names.Estimated{i})
            title(Res.Parameters.Names.Estimated{i})
        end
    end
else
    k = ceil(sqrt(dim));
    for i = 1:dim  
        subplot(k,k,i)
        plot(Res.Thetas(i,:))
        title(Res.Parameters.Names.Estimated{i})
    end
end
    
figure(2)
subplot(2,1,1)
plot(Res.LogLiks)
title('LogLiks')
subplot(2,1,2)
plot(Res.LogPosts)
title('LogPosts')

dim = Res.Parameters.NbParsEstimated;
figure(3)
clf
if dim <6
    for i = 1:dim
        for j = i+1:dim
            subplot(dim,dim,(i-1)*dim+j)
            plot(Res.Thetas(j,:),Res.Thetas(i,:),'.')
            if i == 1
                title(Res.Parameters.Names.Estimated{j})
            end
            subplot(dim,dim,(j-1)*dim+i)
            plot(Res.TransfThetas(i,:),Res.TransfThetas(j,:),'.')
            if i == 1
                ylabel(Res.Parameters.Names.Estimated{j})
            end
        end
        subplot(dim,dim,(i-1)*dim+i)
        plot(Res.TransfThetas(i,:))
        if i == 1
            ylabel(Res.Parameters.Names.Estimated{i})
            title(Res.Parameters.Names.Estimated{i})
        end
    end
else
    k = ceil(sqrt(dim));
    for i = 1:dim  
        subplot(k,k,i)
        plot(Res.TransfThetas(i,:))
        title(Res.Parameters.Names.Estimated{i})
    end
end

figure(4)
k = ceil(sqrt(dim));
for i = 1:dim 
    subplot(k,k,i)
    temp = AutoCorrelation(Res.TransfThetas(i,end-NbSamples+1:end),MaxLag);
    plot(temp)
    ylim([-0.2 1])
    title(length(Res.TransfThetas(i,:))/(1+2*sum(temp(2:end))))
end