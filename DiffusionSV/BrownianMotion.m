function [B]=BrownianMotion(t_1,t_2,m)

% Generates m+1 values of a Brownian motion between 
% points t_1 and t_2 with initial value B(t_1)=0

% NOTE: B(1),...,B(m+1) corresponds to B(t_1),...,B(t_2) resp.

B=zeros(1,m+1);
B(1) = 0;              % B(1) is the initial value of the B-motion
step = (t_2-t_1)/m;  
vol = 1.0;
B(2:m+1)=normrnd(0,sqrt(vol*step),1,m);
B=cumsum(B);