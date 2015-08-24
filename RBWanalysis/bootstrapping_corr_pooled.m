% bootsrapping correlations:

% draw 2 random samples from the pooled data( control and gradient e.g.)and
%correlate, iterate X times

iterations=10000;
ssize=1000;
tic

correlation=NaN;
pvalue=NaN;

for i= 1:iterations
    
    S1=NaN(ssize,1);
    S2=NaN(ssize,1);
    
    for bin= 1:(size(datapool,2)) %draw samples for each bearing bin
        
        S1(:,bin) = randsample(datapool{bin},ssize);
        S2(:,bin) = randsample(datapool{bin},ssize);
        
    end
    
    [correlation(i),pvalue(i)]=corr(nanmean(S1)',(mB)');
    
end
toc
%% plot histogram and variance:

%variance of distribution:
phi=mean(correlation)+(abs(std(correlation))*2);
phiidx=find(round(b*50)/50==round(phi*50)/50);


[v,b]=normhist(correlation,20,0);
hold on
name=dirname2(cd);
title(['pooled: ' name ' # of iterations: ' num2str(iterations)])
plot([phiidx phiidx ],[0 0.15],'--k')
text(phiidx,0.15,['2phi=' num2str(round(phi*1000)/1000)])