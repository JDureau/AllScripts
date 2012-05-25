function [] = PlotRes2(Res,Parameters)

Names = Parameters.EstimatedParameterNames;
k = length(Names);
clf
for i = 1:k
    subplot(k,k,(i-1)*k + i)
    plot(Res.Thetas(i,:))
end
for i = 1:k-1
    for j = i+1:k
        subplot(k,k,(i-1)*k + j)
        scattercloud(Res.Thetas(j,:),Res.Thetas(i,:))
    end
end
for i = 2:k
    for j = 1:i-1
        subplot(k,k,(i-1)*k + j)
        scattercloud(Res.LogThetas(j,:),Res.LogThetas(i,:))
    end
end
for i = 1:k
    subplot(k,k,i)
    title(Names{i})
end
for i = 1:k
    subplot(k,k,(i-1)*k+1)
    ylabel(Names{i})
end
