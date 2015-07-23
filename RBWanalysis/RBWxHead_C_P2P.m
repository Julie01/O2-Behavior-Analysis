
% code bins eigenworm analysis variable by bearing
% for this code all als-files ans EigenWormTracks-files must be in a folder, one for each experiment
% -the GT loop will run over control and gradient folders
% -The sumAll parameter will sum all angles if set to 1, if 0 it will take
% -the segement defined wih parameter segment
% -peakcutoff defines minimum rbw peak amplitude

clear
set(0,'DefaultTextInterpreter','none');
home=cd;
name=dirname(cd);
binningtype= 'RBW';
sumAll=0;
segment=2;
peakcutoff=0.35;

fig=figure;

lo=1;
hi=179;


for GT =1:2 % if you have 2 folders with control and gradient data
    tic
    if GT==1
        cd('control')
        
    elseif GT==2
        cd('gradient')
    end
    
    %%%%
    
    mRBW=NaN(15,20);
    mB=NaN(15,20);
    mC=NaN(15,20);
    
    files=dir('*_als.mat');
    EWfiles=dir('*als_V8_EigenWormTracks.mat');
    
    % create gradient matrix
    gradient=zeros(1060,2060);
    gl=21:-17/(2000):4;
    low=4;
    high=21;
    
    EC=1;
    
    for i=1:size(gradient,1)
        gradient(i,1:30)=high;
        gradient(i,31:2031)=gl;
        gradient(i,2032:end)=low;
    end
    
    
    %% do analysis for various variables and put in subplots
    Vidx= 2  ; % variable index, 2= RBW,
    
%     for E=[ 1 2 3 4 5 6 8 9] 
    for E=1:length(files)
        
        disp(E)
        tic
        disp('...loading')
        load(files(E).name)
        load([files(E).name(1:end-4) '_anglesV7.mat']);
        
        load(EWfiles(E).name)
        toc
        
        FNames = fieldnames(EigenWormTracks);
        
        %     FNames{Vidx}
        
        p_hC=NaN;
        AllPeakHC=[];
        AllPeakCurving=[];
        AllPeakRBW=[];
        %%
        for TN=1:length(Tracks) % go through all tracks of current video
            
            % get head positions:
            BBy=(round(Tracks(1,TN).BoundingBox(2:4:end))*10)/10;
            BBx=(round(Tracks(1,TN).BoundingBox(1:4:end))*10)/10;
            Head=WormGeomTracks(1,TN).WormHeadPosition;
            HeadX=medfilt1(WormGeomTracks(1,TN).WormHeadPosition(:,1)+BBx',4);
            HeadY=medfilt1(WormGeomTracks(1,TN).WormHeadPosition(:,2)+BBy',4);
            
            X=  (Tracks(1,TN).SmoothX);
            Y=  (Tracks(1,TN).SmoothY);
            
            %(1) check if track has unusual high angular speed or low velocity
            
            meanAngVelocity=nanmeanJ(abs(Tracks(1,TN).AngSpeed));
            meanSpeed=nanmeanJ(abs(Tracks(1,TN).Speed));
            
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
                
                [bearing,curving]=getbearing_d15(X,Y,100,0,Tracks,TN);
                bearing=bearing';
                
                bi= find( X<600 | X>1800 | Y<160 | Y>970 | bearing < lo | bearing > hi); %%deletes data which are outside central arena
                nx=[bi,nx];
                
                HeadX(nx)=NaN;
                HeadY(nx)=NaN;
                
                C_ind=sub2ind(size(gradient),HeadY,HeadX);
                
                %get O2 concentration at this position of path:
                nan_idx=isnan(C_ind);
                C_ind=C_ind(nan_idx~=1);
                
                %get O2 concentration at this position of path and reinsert Nans:
                sensoryPathNaN=(gradient(round(C_ind)));
                sensPath=NaN(1,length(nan_idx));
                sensPath(nan_idx~=1)=sensoryPathNaN;
                C=sensPath; % C / sec  for 10hz data
                
                
                %%%
                rbw = EigenWormTracks(TN).(FNames{Vidx});
                if sumAll==0
                    rbw = rbw(:,segment)';
                elseif sumAll==1
                    rbw = sum(rbw,2);
                end
                
                
                rbw(nx)=NaN;
                %%%
                
                %find RBW peaks and get corresponding bearing and
                %curving values:
                
                [maxR minR]=peakdet(rbw,0.40);
                
                
                if ~isempty(maxR) & ~isempty(minR)
                    
                    gi=maxR(:,2)>peakcutoff;
                    maxR=maxR(gi,:);
                    gi=minR(:,2)<-peakcutoff;
                    minR=minR(gi,:);
                    peaks=vertcat(minR, maxR);
                    
                    [maxR2 minR2]=peakdet(rbw,0.12);
                    smallpeaks=vertcat(minR2, maxR2);
                    
                    %%plot RBW and peaks:
                    %    cla
                    %                 if TN<10
                    %                     figure
                    %                     hold on
                    %                     plot (rbw)
                    %                     scatter(smallpeaks(:,1),smallpeaks(:,2),'c')
                    %                     scatter(maxR(:,1),maxR(:,2),'r')
                    %                     scatter(minR(:,1),minR(:,2),'m')
                    %                 end
                    
                    %%plot RBW, C and peaks:
                    %                 if TN<100
                    %                     figure(fig1)
                    %                     hold on
                    %                     scatter(HeadX,HeadY,10,C)
                    %                     hold on
                    %                     scatter(HeadX+10,HeadY,10,rbw)
                    %                     scatter(HeadX(peaks(:,1))+10,HeadY(peaks(:,1)),'m')
                    %                     caxis([-0.25 0.25])
                    % %
                    %                 end
                    
                    
                    %%%%%%%%%%%%%%
                    %find previous peak/valley &
                    % delta C between from last peak to peak
                    
                    if length(peaks)>3
                        p_hC=NaN;
                        pre_peak=NaN;
                        
                        for p=1:length(peaks)
                            dp = (smallpeaks(:,1)-peaks(p,1));
                            [dps,IX]=sort(dp);
                            [pp,pre_p_idx]=max(dps(dps<0));
                            
                            if ~isempty(pp)
                                pre_peak(p)=smallpeaks(IX(pre_p_idx),1);
                                p_hC(p)=C(peaks(p,1))-C(pre_peak(p));
                            else
                                pre_peak(p)=NaN;
                                p_hC(p)=NaN;
                            end
                            
                            %omit peaks with previous NaN stretch
                            try
                                if length(find(isnan(rbw(pre_peak(p):peaks(p,1)))))>8
                                    pre_peak(p)=NaN;
                                    p_hC(p)=NaN;
                                end
                            end
                            
                        end
                        
                        %                      scatter(pre_peak(~isnan(pre_peak)),rbw(pre_peak(~isnan(pre_peak))),'mp')
                        
                        AllPeakHC=horzcat(AllPeakHC,p_hC);
                        AllPeakRBW=horzcat(AllPeakRBW,peaks(:,2)');
                        
                    end
                    
                end
                
            end
            
        end % end Track loop
        
        % NaN inf and outlier values:
        cutoff = (nanmean(AllPeakHC)) + (nanstd(AllPeakHC)*2.5);
        fi=find(~isfinite(AllPeakHC)| AllPeakHC >cutoff | AllPeakHC<-cutoff);
        AllPeakHC(fi)=NaN;
        
        
        
        %%% binning
        if binningtype=='RBW'
            %% bin by rbw
            AllPeakRBW=abs(AllPeakRBW);
            cc=1;
            %bin -->for what?
            a=0.02;
            for i=0.35:a:0.62
                
                bin_idx=AllPeakRBW>i &AllPeakRBW<=i+a;
                pb1=AllPeakHC(bin_idx);
                %         pc1= AllPeakCurving(bin_idx);
                rbw_bin=abs(AllPeakRBW(bin_idx));
                
                mB(EC,cc)=nanmeanJ(pb1);
                %         mC(E,cc)=nanmeanJ(pc1);
                mRBW(EC,cc)=nanmeanJ(rbw_bin);
                
                cc=cc+1;
                
            end
            
        else
            %%%% bin by HC
            AllPeakRBW=abs(AllPeakRBW);
            cc=1;
            %bin -->for what?
            a=0.04;
            for i=-0.2:a:0.2
                
                bin_idx=AllPeakHC>i & AllPeakHC<i+a;
                pb1=AllPeakHC(bin_idx);
                %         pc1= AllPeakCurving(bin_idx);
                rbw_bin=(AllPeakRBW(bin_idx));
                
                mB(E,cc)=nanmeanJ(pb1);
                %         mC(E,cc)=nanmeanJ(pc1);
                mRBW(E,cc)=nanmeanJ(rbw_bin);
                
                cc=cc+1;
            end
            
        end % end if binning type
        
        EC=EC+1;
        
    end % end files loop
    mB=mB(1:EC-1,1:cc-1);
    mRBW=mRBW(1:EC-1,1:cc-1);
    save mRBW mRBW
    save mB mB
    
    
    %% plotting
    
    if binningtype=='RBW'
        
        %% plot:(1) variable binned by RBW
        
        figure(fig)
        set(fig, 'name',name)
        
        sem=nanstd(mB,1)/sqrt(EC-1);
        hold on
        CM=(winter(2)/1.5);
        plot(nanmeanJ(mB,1),'linewidth', 2, 'color' ,CM(GT,:));
        errorb(nanmeanJ(mB,1),sem)
        set(gca,'XTick',1:length(mRBW))
        set(gca,'XTickLabel',(round(nanmeanJ(mRBW)*100)/100));
        xlabel('peak RBW')
        ylabel('Head dc (%O2/s)')
        title (['dC_vs_RBW_head  ' num2str(lo) '-' num2str(hi) ' deg ..'  name])
         ylim([-0.1 0.1])
        ylim auto
        
        
    else
        
        %% plot:(2) variable binned by HC
        
        
        figure(fig)
        set(fig, 'name',name)
        
        sem=nanstd(mRBW,1)/sqrt(EC-1);
        hold on
        CM=(winter(2)/1.5);
        plot(nanmeanJ(mRBW,1),'linewidth', 2, 'color' ,CM(GT,:));
        errorb(nanmeanJ(mRBW,1),sem)
        set(gca,'XTick',1:length(mB))
        set(gca,'XTickLabel',(round(nanmeanJ(mB)*1000)/1000));
        ylabel('peak RBW')
        xlabel('Head dC (%O2/s)')
        title (['C_vs_RBW_head  ' num2str(lo) '-' num2str(hi) ' deg ..'  name])
        
        
        %%%%
    end
    cd(home)
    toc
end % end folder loop

%%

legend('control','gradient ')
if sumAll==1
    saveas(gcf,['dC_vs_RBW_sum ' num2str(lo) '-' num2str(hi) '.fig'])
    disp('sum')
else
    saveas(gcf,['dC_vs_RBW_head ' num2str(lo) '-' num2str(hi) '.fig'])
end
disp('done')
y = chirp( 0:0.0008:2,5);
sound(y);




