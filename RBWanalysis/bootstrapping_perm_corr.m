% bootsrapping correlations:

% draw 2 random samples from the pooled data( control and gradient e.g.)and
%correlate, iterate X times

iterations=1000000;
ar=input('actual r?');

correlation=NaN;
pvalue=NaN;
mB=5 :5:175;

for i= 1:iterations
    
    S1 = randsample(mB,length(mB));% shuffle bearing bins
    
    [correlation(i),pvalue(i)]=corr(nanmean(mRBW)',S1');
    
end

%% plot histogram:
%variance of distribution:
figure
rf=50;
[v,b]=normhist(correlation,20,0);
phi=mean(correlation)+(abs(std(correlation))*2);
phiidx=find(round(b*rf)/rf==round(phi*rf)/rf);
while isempty(phiidx)
    rf=rf*0.9;
    phiidx=find(round(b*rf)/rf==round(phi*rf)/rf);
end
phi2=mean(correlation)-(abs(std(correlation))*2);
phiidx2=find(round(b*rf)/rf==round(phi2*rf)/rf);
while isempty(phiidx2)
    rf=rf*0.9;
    phiidx2=find(round(b*rf)/rf==round(phi2*rf)/rf);
end


hold on
name=dirname2(cd);
name2=dirname(cd);
title([ name '-' name2 ' # of iterations: ' num2str(iterations)])
plot([phiidx phiidx ],[0 0.10],'--k')
plot([phiidx2 phiidx2 ],[0 0.10],'--k')
text(phiidx,0.10,['2phi=' num2str(round(phi*1000)/1000)])


if ar~=0
ridx=find(round(b*rf)/rf==round(ar*rf)/rf);
while isempty(ridx)
    rf=rf*0.9;
    ridx=find(round(b*rf)/rf==round(ar*rf)/rf);
end
l=plot([ridx ridx ],[0 0.12],'r');
text(ridx,0.12,['actual r=' num2str(round(ar*1000)/1000)])
end

%p values:

if ar>mean(correlation)
p_gradient=length(find(correlation>ar))/length(find(correlation<ar));
else
    p_gradient=length(find(correlation<ar))/length(find(correlation>ar));
end

legend(l,['p=' num2str(p_gradient)])

