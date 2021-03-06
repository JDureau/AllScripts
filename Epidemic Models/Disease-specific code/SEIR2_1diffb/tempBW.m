

m = 5;
x0 = 0.2;
asympt = 1;
k = 2;
tinfl =  1/k*log(B*(1/(1-m)));


tis = 0:0.1:600;
x0 = 0.5+0.2*rand(1,1);
A = 0.9;
m = 10.1+15*rand(1,1);
% m = 0.9*rand(1,1);
tinfl = 500 + 100*rand(1,1);

B = (1 - (x0/A)^(1-m));
k = abs(1/tinfl*log(B/(1-m)));

NbParts = 1;
xis = [];
xis(1:NbParts,1) = (x0^(1-m)-A^(1-m))/(1-m);
sigma = 5;
for i = 2:length(tis)
    xis(1:NbParts,i) = xis(1:NbParts,i-1) -k*xis(1:NbParts,i-1)*0.1 + sigma*xis(1:NbParts,i-1)*0.01.*randn(NbParts,1);
end
betas =  ((1-m)*xis+A^(1-m)).^(1/(1-m));
plot(tis,betas')
ylim([0 1])
title(m)

B = (1 - (x0/A)^(1-m));
k = abs(1/tinfl*log(B/(1-m)));
k

clf
xs = A*(1-B*exp(-k*tis)).^(1/(1-m));
plot(tis,xs)
ylim([0 1])
title(m)







% clear x
% tmp = solve('exp(k*tinfl)*(1-x)+(x0/asympt)^(1-x)-1=0')
% 
% tmp = solve('5*(1-x)+4^(1-x)-1=0')
% double(tmp)


a = x0/asympt;
b = exp(k*tinfl);
tmp = 1 - 1/b + lambertw(2,(a^(1 - (-1 + b)/b)* log(a))/b)/log(a);
real(tmp)

A = asympt;
B = (1 - (x0/asympt)^(1-m));

ts = -10:0.01:10;
xs = A*(1-B*exp(-k*ts)).^(1/(1-m));

plot(ts,xs)


[b,ind] = max(diff(xs));
'tinfl' 
ts(ind)

 
1/k*log(B*(1/(1-m)))


tinfl = 400;
k = 1;

A = asympt;

tmp = solve('(30)^(1-x)+390*(1-x)-1=0');
m = double(tmp)

B = (A^(1-m) - x0^(1-m))/(A^(1-m));


ms = 0.01:0.01:2;
ys = (x0/asympt).^(1-ms)+k*tinfl*(1-ms)-1;

plot(ms,ys)



x =  0.01:0.001:20;
y = []
for i = 1:length(x)
    y(i) = (x0/A)^(1-x(i))+exp(k*tinfl)*(1-x(i))-1;
end
plot(x,log(y))



