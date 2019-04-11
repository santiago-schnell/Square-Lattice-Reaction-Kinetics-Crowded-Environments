B=6*ones(N,N,N);

for type = 1:5
    for mol = 1:n(type)
        B(M(mol,type,1),M(mol,type,2),M(mol,type,3))=type;
    end
end

error = 0;
for i=1:N
    error = error + norm(A(:,:,i)-B(:,:,i));
end


