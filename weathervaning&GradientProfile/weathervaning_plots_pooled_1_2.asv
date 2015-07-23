%% weathervaning bar plots
% use on runinfo files
clear

%plot?
plot=1;
files =dir('*runinfo*');

%parameters:
first=5;
last=175;
binsize=20;
numbins=180/binsize;

% for each experiment: put all run data into one vector:
mP= NaN(30,numbins);
mT= NaN(30,numbins);
mS= NaN(30,numbins);
c=0;

for batch=[1 3]
    
    load(files(batch).name);
    disp(files(batch).name);
   

for F=1:length(bearing)
    
    if ~isempty(bearing{F})
    
    c=c+1;
    bearingAll=[];
    dBearingAll=[];
    SpeedAll=[];

    
    for i=1:size(bearing{F},1)
        if ~isempty(dBearing{F}{i,1})
            for rr=1:size(bearing{F}(i,:),2)
                bv=bearing{F}{i,rr};
                tr=dBearing{F}{i,rr};
                sp=speed_ctr{F}{i,rr};
               
                bearingAll=cat(2,bearingAll,bv);
                dBearingAll=cat(2,dBearingAll,tr);
                SpeedAll=cat(2,SpeedAll,sp);

            end
        end
    end
    
    
    bv=round(bearingAll*10)/10;
    tr=round(dBearingAll*10)/10;
    sp=round(SpeedAll*100)/100;

    %kill extremely high turning rates (short reversals etc)
    
    bi=find(tr>19 | tr<-19);
    bi2= find(bv>last| bv<first);
    bi=[bi,bi2];
    tr(bi)=NaN;
    bv(bi)=NaN;
    sp(bi)=NaN;
    
if ~isempty(bv)
    
    [X,v]=hist(bv,10);
    X1(c,:)=X./sum(X);
    cc=1;
    %bin -->for what?
    
    for i=0 :binsize:180
        bin_idx= bv>i & bv<=i+binsize;
        pb1=bv(bin_idx);
        pt=tr(bin_idx);
        st=sp(bin_idx);

        mP(c,cc)=nanmean(pb1);
        mT(c,cc)=nanmean(pt);
        mS(c,cc)=nanmean(st);

        cc=cc+1;
    end
    
end
    end
end

end % end files loop


nd=(cd);
d= strfind(cd, '\');
name=nd(d(end-1)+1:end);


%% plot:(1) turning rate binned by bearing
if plot==1

figure
sem=nanstd(mT*3,1)/sqrt(c);
hold on
bar(nanmean(mT*3,1));
errorb(nanmean(mT*3,1),sem)
set(gca,'XTick',[1:length(mP)])
set(gca,'XTickLabel',round(nanmedian(mP,1)));
xlabel('bearing')
ylabel('turning bias (deg/sec)')
ylim([-1 0.5])
title (name)

if ~exist ('plots', 'dir')
mkdir('plots')
end
cd('plots')
saveas(gca, 'weathervaning_pooled_olddata.fig')
% saveas(gca, 'weathervaning.jpg')
cd ..\

save turningbias mT



end


%% (2) bearing  binned by turning rate


