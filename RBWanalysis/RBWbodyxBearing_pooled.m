
% bins eigenworm analysis variable by bearing

clear
home=cd;
fig=figure;
segment=11;
peakcutoff=0.075;
peakdelta=0.15;
meanRBW=NaN(2,35);

for GT =1:2 % if you have 2 folders with control and gradient data
    
    if GT==1
        cd('control')
        
    elseif GT==2
        cd('gradient')
    end
     
    AllPeakBearing=[];
    AllPeakCurving=[];
    AllPeakRBW=[];
    
    files=dir('*_als.mat');
    EWfiles=dir('*als_V8_EigenWormTracks.mat');
    
    %% do analysis for various variables and put in subplots
    Vidx=  [ 2 ] ; % variable index, 2= RBW,
    
%     for E=1:11
        for E=1:length(files)
        E
        disp('...loading')
        load(files(E).name)
        load(EWfiles(E).name)
        FNames = fieldnames(EigenWormTracks);
        
        
        %%
        for TN=1:length(Tracks) % go through all tracks of current video
            
            RunX=  (Tracks(1,TN).SmoothX);
            RunY=  (Tracks(1,TN).SmoothY);
            
            %(1) check if track has unusual high angular speed or low velocity
            
            meanAngVelocity=nanmean(abs(Tracks(1,TN).AngSpeed));
            meanSpeed=nanmean(abs(Tracks(1,TN).Speed));
            
            if meanAngVelocity<30 & meanSpeed>0.04 % continue only if speed is not high and angular speed is low
                
                %(2) remove omegas and reversals:
                Rev1=[];
                Omegas=[];
                
                for ii =1:size(Tracks(1,TN).polishedReversals,1)
                    Rev=(Tracks(1,TN).polishedReversals(ii,1):Tracks(1,TN).polishedReversals(ii,2));
                    Rev1=horzcat(Rev1,Rev);
                end
                
                for ii =1:size(Tracks(1,TN).OmegaTrans,1)
                    Omegas1=(Tracks(1,TN).OmegaTrans(ii,1):Tracks(1,TN).OmegaTrans(ii,2));
                    Omegas=horzcat(Omegas,Omegas1);
                end
                nx=horzcat(Rev1,Omegas);
                
                [bearing,curving]=getbearing_d15(RunX,RunY,100,0,Tracks,TN);
                
                bi= find( RunX<600 | RunX>1800 | RunY<160 | RunY>970); %%deletes data which are outside central arena
                nx=[bi,nx];
                bearing(nx)=NaN;
                curving(nx)=NaN;
                %%%
                ctrl=abs(EigenWormTracks(TN).RBW(:,11))';
                
                if Vidx<4
                    rbw = EigenWormTracks(TN).(FNames{Vidx});
                    rbw = rbw(:,segment)';
                else
                    rbw = EigenWormTracks(TN).(FNames{Vidx})';
                end
                
                rbw(nx)=NaN;
                %%%
                
                %find RBW peaks and get corresponding bearing and
                %curving values:
                
                [maxR minR]=peakdet(rbw,peakdelta);
                
                
                if ~isempty(maxR) & ~isempty(minR)
                    
                    gi=maxR(:,2)>peakcutoff;
                    maxR=maxR(gi,:);
                    gi=minR(:,2)<-peakcutoff;
                    minR=minR(gi,:);
                    peaks=vertcat(minR, maxR);
                    
                    %%plot RBW and peaks:
                    if GT==1 & E==1 & TN<1
                        figure
                        cla
                        hold on
                        plot (rbw)
                        scatter(maxR(:,1),maxR(:,2),'r')
                        scatter(minR(:,1),minR(:,2),'c')
                        title (['RBW body peak cutoff/delta=' num2str(peakcutoff) '/' num2str(peakdelta)])
                    end
                    
                    if size(bearing,2)==1
                        bearing=bearing';
                    end
                    
                    p_bearing=bearing(peaks(:,1));
                    p_curving=curving(peaks(:,1));
                    
                    
                    AllPeakBearing=horzcat(AllPeakBearing,p_bearing);
                    AllPeakRBW=horzcat(AllPeakRBW,peaks(:,2)');
                    AllPeakCurving=horzcat(AllPeakCurving,p_curving');
                    
                end
                
            end
            
        end % end Track loop
        
    end % end files loop
    
    %% bin by bearing
    mRBW=NaN(1,4);
    mB=NaN(1,4);
    mC=NaN(1,4);
    cc=1;
    %bin -->for what?
    a=5;
    for i=5 :a:175
        
        bin_idx= AllPeakBearing>i & AllPeakBearing<i+a;
        pb1= AllPeakBearing(bin_idx);
        pc1= AllPeakCurving(bin_idx);
        rbw_bin=abs(AllPeakRBW(bin_idx));
        
        mB(1,cc)=nanmean(pb1);
        mC(1,cc)=nanmean(pc1);
        mRBW(1,cc)=nanmean(rbw_bin);
        rbw_bins{cc}=rbw_bin;
        sem=nanstd(rbw_bin,1)/sqrt(E);
        
        cc=cc+1;
    end
    
    save rbw_bins rbw_bins
    

    %% plot:(1) variable binned by bearing
    name=dirname(cd);
    CM=(winter(2)/1.5);
    figure(fig)
    set(fig, 'name',name)
    
    hold on
    scatter(1:cc-1,mRBW,'markerfacecolor',CM(GT,:));
    %errorb(nanmean(mRBW,1),sem,'color',[0.5 0.5 0.5])
    set(gca,'XTick',[1:2:length(mB)])
    set(gca,'XTickLabel',(round(mB(1:2:end)*1)/1));
    xlabel('bearing')
    ylabel('mean rbw peak value')
    title (['RBW body peak cutoff/delta=' num2str(peakcutoff) '/' num2str(peakdelta)])
    ylim auto
    ylim([0.12 0.139])
    
    [h(GT),p(GT)]=corr(mRBW',mB')
    
    
    %%%%
    
    cd(home)
    
end % end folder loop
%%
legend(['control r=' num2str([h(1),p(1)])],['gradient r=' num2str([h(2),p(2)])])

disp('done')




