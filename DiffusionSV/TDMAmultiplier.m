function d = TDMAmultiplier(a,b,c,x)
%a, b, c are the column vectors for the compressed tridiagonal matrix A
%x is the vector, d=A*X
%Note by convention a(1)=c(n)=0
n = length(b); % n is the number of rows
d=b.*x;
d(1)=d(1)+c(1)*x(2);
d(n)=d(n)+a(n)*x(n-1);
d(2:n-1)=d(2:end-1)+a(2:n-1).*x(1:n-2)+c(2:n-1).*x(3:n);
end
