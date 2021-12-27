% load test data: urban data set, initial values and test parameters
load data.mat

% create a first order roughness penalty matrix and add dicontinuities at
% bands 100 and 120
D=createD(M);
dc=[100 120];
D(dc,:)=0;



[A,S,Jout]=unmixing_lq(q,Y,A0,S0,h,gam,D,del,1000);

for i=1:r
    subplot(r,2,2*i-1);
    plot(A(:,i));
    title(['endmember ' num2str(i)])
    subplot(r,2,2*i);
    imagesc(reshape(S(i,:),307,307))
    title(['abundance map ' num2str(i)])

end    
