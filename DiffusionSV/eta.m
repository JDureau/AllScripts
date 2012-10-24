function [Xd]=eta(BB,dtj,m,Yd1,Yd2)

kk=0:m+1;
s = dtj*kk/(m+1);

Xd = s * Yd2 + ( dtj - s ) * Yd1;

Xd = BB + Xd/dtj;	

