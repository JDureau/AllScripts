function xout = logit(xin)

if or(xin<=0,xin>=1)
    disp('pb')
end

xout = log(xin./(1-xin));

