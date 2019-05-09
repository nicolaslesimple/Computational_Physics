function val = off(A)

N=length(A);
val = 0;
for i=1:N
    for j=1:N
        if (i~=j)
            val=A(i,j)*A(i,j)+val;
        end
    end
end
val=sqrt(val);
end