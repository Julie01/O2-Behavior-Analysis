% bootsrapping correlations:

% draw 2 random samples from the pooled data( control and gradient e.g.)and
%correlate, iterate X times

iterations=1000000;
ssize=11;
correlation=NaN;
pvalue=NaN;

cr=input('control r?');
ar=input('gradient r?');

tic

%load data
cd('control')
load('mRBWb');
d1=mRBW;
cd ..\
cd('gradient')
load('mRBWb');
cd ..\
d2=mRBW;
datapool=vertcat(d1,d2);

mB=5 :5:175;


for i= 1:iterations
    
    S1=NaN(11,1);
    
    for bin= 1:(size(datapool,2)) %draw samples for each bearing bin
        
        S1(:,bin) = randsample(datapool(:,bin),ssize,1);
        
    end
    
    [correlation(i),pvalue(i)]=corr(nanmean(S1)',(mB)');
    
end

toc
%% plot histogram:
%variance of distribution:
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
title([ name ' # of iterations: ' num2str(iterations)])
plot([phiidx phiidx ],[0 0.12],'--k')
plot([phiidx2 phiidx2 ],[0 0.12],'--k')
text(phiidx,0.10,['2phi=' num2str(round(phi*1000)/1000)])


if ar~=0
ridx=find(round(b*rf)/rf==round(ar*rf)/rf);
while isempty(ridx)
    rf=rf*0.9;
    ridx=find(round(b*rf)/rf==round(ar*rf)/rf);
end
l(1)=plot([ridx ridx ],[0 0.12],'r');
text(ridx,0.12,['gradient r=' num2str(round(ar*1000)/1000)])
end


if cr~=0
ridx=find(round(b*rf)/rf==round(cr*rf)/rf);
while isempty(ridx)
    rf=rf*0.9;
    ridx=find(round(b*rf)/rf==round(cr*rf)/rf);
end
l(2)=plot([ridx ridx ],[0 0.12],'color', [0.5 0.5 0.5]);
text(ridx,0.12,['control r=' num2str(round(cr*1000)/1000)])
end

%p values:

if ar>mean(correlation)
p_gradient=length(find(correlation>ar))/length(find(correlation<ar));
else
    p_gradient=length(find(correlation<ar))/length(find(correlation>ar));
end
    
if cr>mean(correlation)
p_control=length(find(correlation>cr))/length(find(correlation<cr));
else
    p_control=length(find(correlation<cr))/length(find(correlation>cr));
end

legend(l,['p=' num2str(p_gradient)],['p=' num2str(p_control)])
