
% bins eigenworm analysis variable by bearing

clear
home=cd;
fig=figure;
segment=2;
peakcutoff=0.15;
peakdelta=0.2;
meanRBW=NaN(2,35);

for GT =1:2 % if you have 2 folders with control and gradient data
    
    if GT==1
        cd('control')
    elseif GT==2
        cd('gradient')
    end
    
    mRBW=NaN(2,4);
    mB=NaN(2,4);
    mC=NaN(2,4);
    
    files=dir('*_als.mat');
    EWfiles=dir('*als_V8_EigenWormTracks.mat');
    
    %% do analysis for various variables and put in subplots
    Vidx=  [ 2 ] ; % variable index, 2= RBW,
    
    % for E=[ 1 2 3 4 5 6 8 9]
    for E=1:length(files)
        E
        disp('...loading')
        load(files(E).name)
        load(EWfiles(E).name)
        
        FNames = fieldnames(EigenWormTracks);
        
        FNames{Vidx}
        
        AllPeakBearing=[];
        AllPeakCurving=[];
        AllPeakRBW=[];
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
                    if E==1 & GT==1 & TN< 10 & length(rbw)>3000
                    figure
                    cla
                    hold on
                    plot (rbw)
                    scatter(maxR(:,1),maxR(:,2),'r')
                    scatter(minR(:,1),minR(:,2),'c')
                    title ([num2str(TN) ': RBW head peak cutoff/delta=' num2str(peakcutoff) '/' num2str(peakdelta)])
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
        
        
        
        %% bin by bearing
        cc=1;
        %bin -->for what?
        a=5;
        for i=5 :a:175
            
            bin_idx= AllPeakBearing>i & AllPeakBearing<i+a;
            pb1= AllPeakBearing(bin_idx);
            pc1= AllPeakCurving(bin_idx);
            rbw_bin=abs(AllPeakRBW(bin_idx));
            
            mB(E,cc)=nanmean(pb1);
            mC(E,cc)=nanmean(pc1);
            mRBW(E,cc)=nanmean(rbw_bin);
            
            cc=cc+1;
        end
        save mRBW mRBW
        
    end % end files loop
    
    meanRBW(GT,:)=nanmean(mRBW,1);
    
    %% plot:(1) variable binned by bearing
    name=dirname(cd);
    
    figure(fig)
    set(fig, 'name',name)
    
    sem=nanstd(mRBW,1)/sqrt(E);
    hold on
    CM=(winter(2)/1.5);
    scatter(1:cc-1,nanmean(mRBW,1),'markerfacecolor',CM(GT,:));
    set(gca,'XTick',[1:length(mB)])
    set(gca,'XTickLabel',(round(nanmedian(mB)*1)/1));
    xlabel('bearing')
    ylabel('mean rbw peak value')
    title (['RBW head peak cutoff/delta=' num2str(peakcutoff) '/' num2str(peakdelta)])
    ylim auto
%    ylim([peakcutoff+0.05 peakcutoff+0.1])
    
    [h(GT),p(GT)]=corr(nanmean(mRBW)',nanmean(mB)');
    %%%%
    
    cd(home)
    
end % end folder loop

legend(['control r=' num2str([h(1),p(1)])],['gradient r=' num2str([h(2),p(2)])])

disp('done')




