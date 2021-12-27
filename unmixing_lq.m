function [A,S,Jout]=unmixing_lq(q,Y,A0,S0,h,gam,D,del,numIters)
% Usage: [A,S,Jout]=unmixing_lq(q,Y,A0,S0,h,gam,D,del,numIters)
% 
% Input: 
%         q: The norm 0<=q<=1
%         Y: The hyperspectral data matrix
%        A0: Initial value for the endmember matrix
%        S0: Initial value for the abundances matrix
%        h:  sparseness tuning parameter
%       gam: smoothness tuning parameter
%         D: the smoothness matrix
%       del: sum to one parameter
%      
% Output:
%        
%         A: The endmember matrix
%         S: The abundances matrix
%         Jout: cost function
% 
% 
%
% J. Sigurdsson, M. O. Ulfarsson and J. R. Sveinsson, 
% "Hyperspectral Unmixing With lq Regularization," in 
% IEEE Transactions on Geoscience and Remote Sensing, vol. 52, 
% no. 11, pp. 6793-6806, Nov. 2014, doi: 10.1109/TGRS.2014.2303155.

    if q>0
        h=h^(2-q);
    else
        h=h^2/2;
    end


    tol=0.00000001;
    maxInnerIt=5000;
    maxIter=10000;
    iterCheck1=5;
    iterCheck2=5;
    
    R=D'*D; 
    gRp=sparse(R.*(R>0)*gam);
    gRn=sparse(-R.*(R<0)*gam);
    Jout=[];
    [M,~]=size(A0);
    [r,P]=size(S0);
    J0=cost(q,Y,A0,S0,h,gam,del,R);
    Yf=[Y; del*ones(1,P)];
    S=S0;Sold=S;
    A1=A0;
    if exist('numIters','var')
        if numIters>0
            tol=-1;
            maxIter=numIters;
        end
    end
    for i=1:maxIter
        i
        Af=[A0; del*ones(1,r)];
        AftAf=Af'*Af;
        try
            evals=eig(AftAf);
        catch
             disp('error evaluating eig(ATA)');
             S=zeros(size(S));
             A=zeros(size(A1));
             return;
        end
        al=max(evals);
        if(al==0)
            error('max eigenvalue of A is zero'); 
        end
        
        lambda=h/(2*al);
        Ba2=0;
        ha=0;
        if lambda>0
            Ba=(2*lambda*(1-q))^(1/(2-q));
            Ba2=lambda*q*exp((q-1)*log(Ba));
            ha=Ba+lambda*q*Ba^(q-1);
        end
        if q==1
            Ba2=lambda;
            ha=lambda;
        end
        alAtA=eye(r)-1/al*AftAf;
        AfTYf=1/al*Af'*Yf;
        Js=J0;Jsiter=1;
        for iter1=1:maxInnerIt
            if q==1 % L1
                S=soft(S+AfTYf-1/al*AftAf*S,lambda);
                S=S.*(S>0);
            else
                Z=AfTYf+alAtA*S;
                if q==0 % L0
                    S=Z.*(Z>ha);
                else % Lq
                    Z=AfTYf+alAtA*S;
                    absZ=abs(Z);
                    % Goran & Solo say 2 iterations is enough
                    Bstar=max(absZ-Ba2,10^-12);
                    for k=1:2
                        Bstar=absZ-lambda*q*exp((q-1)*log(Bstar));
                    end
                    S=sign(Z).*Bstar;
                    S=S.*(absZ>=ha);
                    S=S.*(S>0);
                end
            end
            if rem(iter1,iterCheck1)==0
                Jsiter=Jsiter+1;
                Js(Jsiter)=cost(q,Y,A1,S,h,gam,del,R);
                di=Js(Jsiter-1)-Js(Jsiter);
                di=di/Js(Jsiter);
                if abs(di)<1e-3
                    break;
                end
            end
        end
        if gam>0
            SST=S*S';
            b=-Y*S';
            bp2=b.^2;
            for iter2=1:maxInnerIt
                Cpv=2*(A1*SST+gRp*A1);
                Cnv=2*(gRn*A1);
                Q=(-b+sqrt(bp2+Cpv.*Cnv))./(Cpv);
                Q(isinf(Q))=1;
                A1=A1.*Q;
                if rem(iter2,iterCheck2)==0
                    Mq=mean(abs(Q(A1>10^-5)-1));
                    if Mq<1e-6
                        break;
                    end
                end
            end
        end
        J1=cost(q,Y,A1,S,h,gam,del,R);
        A0=A1;
        J0=J1;
        Jout(i)=J1;
        Ss(i)=sum(abs(S(:)-Sold(:)))/sum(abs(S(:)));     
        if Ss(i)<tol
            break;
        end
        Aold=A1;
        Sold=S;
    end
    A=A1;
end

function [J]=cost(q,Y,A,S,h,gam,del,R)
    [r,P]=size(S);
    AtDtDA=A'*R*A;
    J=0.5*sum(sum((Y-A*S).^2))+0.5*gam*trace(AtDtDA) + ...
        0.5*del^2*sum((ones(1,r)*S-ones(1,P)).^2);
    if q==0
    J=J+0.5*h*sum(S(:)>0);
    else
        %J=J+0.5*h*sum(S(:).^q);%
        J=J+0.5*h*sum(exp(q*log(S(:))));
    end
end

function x_th=soft(x,T)
x_th=sign(x).*max(0,abs(x)-T);
end








