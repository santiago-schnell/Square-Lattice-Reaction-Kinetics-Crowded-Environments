function [M,n]=addmol(M,type,x,y,z,n)

n = n+1;
M(n,type,1) = x;
M(n,type,2) = y;
M(n,type,3) = z;