function xout = invlogit(xin)


xout = exp(xin)./(1+exp(xin));