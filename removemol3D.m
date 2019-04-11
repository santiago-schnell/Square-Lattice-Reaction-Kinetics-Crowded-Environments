% function removes molecule of specified type from specified row of M
% function returns new M and the number of molcules of type n remaining

function [M,n] = removemol(M,type,row,n)

if n < 1, error = 3; end
if type < 1, error = 4; end
if row < 1, error = 5; end

if row < n, M(row:n-1,type,1:3) = M(row+1:n,type,1:3); end

M(n,type,1:3)=0;
n = n-1;

% if norm(M(n+1,:),1) == 0, M = M(1:n,1:10); end
    
