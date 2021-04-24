
function [X]=initialization(N,dim)

    X=rand(N,dim);
    p=rand();
    for i=1:N
        for j=1:dim
            if X(i,j)>p
                X(i,j)=1;
            else
                X(i,j)=0;
            end
        end
    end    
end