function x = TDMAsolver(a,b,c,d)
%a, b, c, and d are the column vectors for the compressed tridiagonal matrix
n = length(b); % n is the number of rows
 
% Modify the first-row coefficients
c(1) = c(1) / b(1);	% Division by zero risk.
d(1) = d(1) / b(1);	% Division by zero would imply a singular matrix. 
 
for i = 2:n
    id = 1 / (b(i) - c(i-1) * a(i));  % Division by zero risk. 
    c(i) = c(i)* id;                % Last value calculated is redundant.
    d(i) = (d(i) - d(i-1) * a(i)) * id;
end
 
% Now back substitute.
x(n) = d(n);
for i = n-1:-1:1
    x(i) = d(i) - c(i) * x(i + 1);
end
end
