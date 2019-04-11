B=6*ones(N,N);

for type = 1:5
    for mol = 1:n(type)
        B(M(mol,type,1),M(mol,type,2))=type;
    end
end

error = norm(A-B);


