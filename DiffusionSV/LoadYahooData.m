function Data = LoadYahooData(file,step)


    tab = csvread(file,1,1);

    nobs = size(tab,1);
    Y = tab(:,1);
    N = (nobs-1)/step;
    
    nobs = length(Y);
    obsstep = N/(nobs-1);
    npoints = N/(nobs-1);
    
    Data = struct();
    Data.Y = Y;
    Data.N = N;
    Data.step = step;
    Data.npoints = npoints;
    Data.nobs = nobs;
    Data.obsstep = obsstep;
    Data.step = step;
    
    
    
%     save([SavePath data_file],'Data');
