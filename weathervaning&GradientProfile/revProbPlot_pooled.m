%% weathervaning hist
% use on runinfo files
% figure
clear
warning off;

CM= jet(6);

for delay=2:4

%plot?
plotting=1;
files =dir('*revinfo*');

%parameters:
first=0;
last=180;
binsize=20;

% for each experiment: put all run data into one vector:
mB= NaN(35,20);
mC= NaN(35,20);
mT= NaN(35,20);
revN= NaN(35,20);
revNnorm=NaN(35,20);
c=0;

for batch=1:length(files)
  
    
    load(files(batch).name);
    disp(files(batch).name);
  
    
   

for F=1:length(bearing)
    
    if ~isempty(bearing{F})
    
    c=c+1;
    bearingAll=[];
    dCdtAll=[];
    SpeedAll=[];
    revAll=[];

    
    for i=1:length(bearing{F})
        if ~isempty(dBearing{F}{i})            
                bv=bearing{F}{i};
%                 dC=diff(sensPath{F}{i});
%                 dC=[dC NaN];
                sp=speed_ctr{F}{i};
                rev=reversals{F}{i};
               
                bearingAll=cat(2,bearingAll,bv);
%                 dCdtAll=cat(2,dCdtAll,dC);
                SpeedAll=cat(2,SpeedAll,sp); 
                revAll=cat(2,revAll,rev); 
        end
    end
    
    
    bv=round(bearingAll*10)/10;
    tr=round(dCdtAll*10)/10;
    sp=round(SpeedAll*100)/100;

    %kill extremely high turning rates (short reversals etc)
    
    bi=find(tr>19 | tr<-19 );
    tr(bi)=NaN;
    bv(bi)=NaN;
    sp(bi)=NaN;
    
    
if ~isempty(bv)
    
    [X,v]=hist(bv,10);
    X1(c,:)=X./sum(X);
    cc=1;
%     delay=15;
    %bin 
    
    for i=first :binsize:(last-binsize)
        
        bin_idx= find(bv>i & bv<=i+binsize);
        bin_idx(bin_idx>length(bv)-delay)=[];
    
        if ~isempty(bin_idx)
%         dC=dCdtAll(bin_idx);
        revV=revAll(bin_idx+delay);
        binB=bv(bin_idx);
        
        mB(c,cc)=nanmean(binB);
%         mC(c,cc)=nanmean(dC);
        revN(c,cc)=nansum(revV);
        revNnorm(c,cc)=nansum(revV)/length(bin_idx);
        end
        
        cc=cc+1;
                
    end
    
end
    end
end

end % end files loop

revNnorm=revNnorm(1:c,1:cc);
% mC=mC(1:c,1:cc);
revN=revN(1:c,1:cc);
mB=mB(1:c,1:cc);

%% plot:(1) turning rate binned by bearing
nd=(cd);
d= strfind(cd, '\');
name=nd(d(end-1):end);

if plotting==1
    
sem=nanstd(revNnorm*3,1)/sqrt(c);
hold on
%plot(nanmean(revNnorm*3,1),'color',CM(delay,:))
shadedErrorBar2([],nanmean(revNnorm*3,1),sem,{'color',CM(delay,:)},0)
set(gca,'XTick',[1:length(mB)])
set(gca,'XTickLabel',round(nanmedian(mB,1)*1)/1);
xlabel('dC/dT')
ylabel('reversal frequency (rev/s)')
ylim([0.005 0.035])
title ([name ])

end

end

saveas(gca, 'reversal_freq.fig')
save revNnorm revNnorm

% save turningbias mT
y=chirp(1:0.001:1.5,30);
sound(y)



