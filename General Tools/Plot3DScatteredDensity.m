function [] = Plot3DScatteredDensity(Data)

if not(size(Data,2) == 2)
    'dimension problem'
    die
end

xmin = min(Data(:,1));
xmax = max(Data(:,1));
ymin = min(Data(:,2));
ymax = max(Data(:,2));

NbSteps = ceil(sqrt(size(Data,1)/8));



XStep = (xmax - xmin)/NbSteps;
YStep = (ymax - ymin)/NbSteps;

res = [];
[X,Y] = meshgrid(xmin+XStep/2:XStep:xmax-XStep/2,ymin+YStep/2:YStep:ymax-YStep/2);
for i = 1:NbSteps
    for j = 1:NbSteps
        indsX = (Data(:,1)>xmin + (i-1)*XStep).*(Data(:,1)<xmin + (i)*XStep);
        indsY = (Data(:,2)>ymin + (j-1)*YStep).*(Data(:,2)<ymin + (j)*YStep);
        inds = indsX.*indsY;
        if sum(inds)
            res(i,j) = length((find(inds)));
        else
            res(i,j) = 0;
        end
    end
end
colormap('default')
contour(X,Y,res')
