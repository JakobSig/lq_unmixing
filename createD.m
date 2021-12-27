function D=createD(T,r)
T=T-1;
D=eye(T);
D=-[D zeros(T,1)]+[zeros(T,1) D(:,1:T)];
if nargin==2
    if r>1
        for i=1:r-1
            D=D'*D;
        end
    end
end
end