% function finds molecule of specified type at specified coords in array
% M and returns its row no. '0' indicates 'molecule not found'.

function row=findmol(M,type,x,y,n)

for i = 1:n
    if M(i,type,1) == x
        if M(i,type,2) == y
            break
        end
    end
end
if i == n, row = 0; end     % FAILURE FLAG - molecule not found.
row=i;


