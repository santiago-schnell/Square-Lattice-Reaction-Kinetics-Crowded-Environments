function [M,n]=addmol(M,type,x,y,n)

n = n+1;
M(n,type,1) = x;
M(n,type,2) = y;