function Data = LoadYahooData(file,step)


    tab = csvread(file,1,1);

    nobs = size(tab,1);
    Y = log(tab(:,1));
%     Y = Y-Y(1);
    N = (nobs-1)/step/253;
    
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
    Data.Z = Sample_Z(N);
    
    
%     save([SavePath data_file],'Data');
