function  KeptSamples = Resample(Weigths)

% [ indx ] = resampleSystematic( Weigths );
% KeptSamples = Samples(indx);
% 
   n = length(Weigths);
   RandU1 = rand(1,1)/n;
   KeptSamples = zeros(1,n);
   CumWeigths = cumsum(Weigths)/sum(Weigths); 
   
   Current = 1;
   for IndResampling = 1:n
       while CumWeigths(Current)<=RandU1+(IndResampling-1)/n
           Current = Current+1;
       end
       KeptSamples(IndResampling) = Current;
   end