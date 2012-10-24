function [BB]=eta_inv(Xd,dt,m)

Yd1=Xd(1);
Yd2=Xd(end);

kk=0:m+1;
s = (kk/(m+1))*dt;

BB = s * Yd2 + ( dt - s ) * Yd1;

BB = Xd - BB/dt;	

