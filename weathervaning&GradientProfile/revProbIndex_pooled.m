%% weathervaning hist
% use on runinfo files

clear
warning off;
home=cd;

for delay=2:4

%plot?
plotting=1;
if plotting==1
    fig=figure;
end

RI=NaN(55,2);

for GT=1:2
    
    if GT==1
        cd('control')
    else
        cd('gradient')
    end
    
    
files =dir('*revinfo*');

%parameters:
first=0;
last=180;
binsize=90;

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

    cc=1;
    
  %bin:  
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

RI(1:size(revNnorm,1),GT)=-1+(revNnorm(:,2)./revNnorm(:,1));
n(GT)=c;
cd(home)

end

%% plot:(1) turning rate binned by bearing

name=dirname2(cd);
cd(home)
%% boxplots:
if plotting==1
    
    figure(fig)
    
    boxplot(RI,'whisker',0.9529,'labels',{'ctrl','gradient'})
    [h,p]=ranksum(RI(1:n(1),1),RI(1:n(2),2));
    text(1.5,0,num2str(h))
    ylim([-0.7 0.75])
    title(name)
    ylabel('relative rev freq change')
    cd('plots')
    saveas(gca, 'reversal modulation boxplot.fig')
    
    cd(home)
    
end

end


saveas(gca, 'reversal_freq.fig')

% save turningbias mT
y=chirp(1:0.001:1.5,30);
sound(y)



