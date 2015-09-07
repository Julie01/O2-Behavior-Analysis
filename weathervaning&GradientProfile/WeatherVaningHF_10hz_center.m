
clear 
% global rois
% cd('mat')
home=cd;
nd=(cd);
d= strfind(cd, '\');
name=nd(d(2):d(end));

files=dir('*als.mat');
nd=(cd);
d= strfind(cd, '\');
name=nd(d(2):d(end));

skelfiles=dir('*angles*');
 % conversion factor for head speed (mm/px) (same as in hr track analyzer)

% cd ..\
for F=1:length(files)
    
    F
    clearvars -except F exp  home name files bearing dBearing speed_ctr...
    HeadSens HeadSpeed HeadSpeedNet
   cc=0;
   cla
    fname=(files(F,1).name)
    load (files(F).name);
    %load([fname(1:end-4), '_anglesV8.mat']);
    
    
    pixToMM=0.0155;
    
    NoRun={};
    Ctr_dBearing=cell(1);
    Ctr_Bearing=cell(1);
    Ctr_Speed=cell(1);
    dSensCtr={NaN};
    dSensHead={NaN};
    HeadRunsX={NaN};
    HeadRunsY={NaN};
    SensoryRunsHead={NaN};
    head_speed={NaN};
    head_speed_net={NaN};
    
    %% ---analyze Tracks----
    
    for T= 1:length(Tracks)
        % only tracks which start after establishment of the gradient
        % (ca frame 500)
        if  Tracks(1,T).Frames(end)>3000 && Tracks(1,T).Frames(1)>1000 
            
            % --plot tracks--
            
            %             hold on
            %             plot( Tracks(1,T).SmoothX, Tracks(1,T).SmoothY,'Color',[0.5 0.5 0.5]);
            %             title(T);
            
           
            
            %% get indices of all behavior events for excluding them from the
            % trajectories:
            
            Rev1=[];
            Omegas=[];
            try
                for ii =1:size(Tracks(1,T).polishedReversals,1)
                    Rev=(Tracks(1,T).polishedReversals(ii,1):Tracks(1,T).polishedReversals(ii,2));
                    Rev1=horzcat(Rev1,Rev);
                end
            catch;end
            
            try
                for ii =1:size(Tracks(1,T).OmegaTrans,1)
                    Omegas1=(Tracks(1,T).OmegaTrans(ii,1):Tracks(1,T).OmegaTrans(ii,2));
                    Omegas=horzcat(Omegas,Omegas1);
                end
            catch;end
            
            NoRun{T}=horzcat(Rev1,Omegas);
            
            RunX=  (Tracks(1,T).SmoothX);
            RunY=  (Tracks(1,T).SmoothY);
            AngVelocity=abs(Tracks(1,T).AngSpeed);
            nx=(NoRun{1,T});
            RunX(nx)=NaN;
            RunY(nx)=NaN;
            
            
            %----exclude those which occur close to border:----
            bi1= find( RunX<600 | RunX>1800);
            bi2= find( RunY<50 | RunY>980);
            bi= setdiff(bi1, bi2);
            bi=[bi,bi2];
            RunX(bi)=NaN;
            RunY(bi)=NaN;
            
            % find start and end points, avoid first and last 4 points at end of run,
            % they often belong already to reversal or omega:
            ci = isnan(RunX);
            RunEnd=find(diff(ci)>0)-5;
            RunStart=find(diff(ci)<0)+5;
            if ci(1)==0;
                RunStart=cat(2,1,RunStart);
            end
            if ci(end)==0;
                RunEnd=cat(2,RunEnd,length(RunX));
            end
            
            
            % find runs longer than X
            X=12;
            longRuns=find((RunEnd-RunStart)>X);
            RunStart=RunStart(longRuns);
            RunEnd=RunEnd(longRuns);
            
            %kill runs which have unusual high angular speed
            %and don't cover 1 small worm travel distance (20 px):
            meanAngVelocity1=NaN(1,1);
            for ii= 1:length(RunStart)
                meanAngVelocity1(ii)=(nanmean(abs(Tracks(1,T).AngSpeed(RunStart(ii):RunEnd(ii)-1))));
                if meanAngVelocity1(ii)>30
                    RunStart(ii)=NaN;
                    RunEnd(ii)=NaN;
                else
                    X=RunX(RunStart(ii):RunEnd(ii));
                    Y=RunY(RunStart(ii):RunEnd(ii));
                    d=NaN(1,1);
                    ci=1;
                    for iii=1:3:length(Y)
                        d(ci)=sqrt(((X(1)-X(iii)).^2)+((Y(1)-Y(iii)).^2));
                        ci=ci+1;
                    end
                    if max(d)<20
                        RunStart(ii)=NaN;
                        RunEnd(ii)=NaN;
                    end
                end
            end
            
            RunStart=RunStart(~isnan(RunStart));
            RunEnd=RunEnd(~isnan(RunEnd));
            RunX1={NaN};
            RunY1={NaN};
            RunFrames1={NaN};
            for ii=1:length(RunStart)
                RunX1{ii}=RunX(RunStart(ii):3:RunEnd(ii)); %downsample
                RunY1{ii}=RunY(RunStart(ii):3:RunEnd(ii));
                RunFrames1{T,ii}=Tracks(1,T).Frames((RunStart(ii):3:RunEnd(ii)));
            end

            %--plot--
            %             if length(RunStart)>1
            %                 for ii= 1:length(RunStart)
            %                     hold on
            %                     %plot([RunX(RunStart(ii)) RunX(RunStart(ii)+6)],[RunY(RunStart(ii)) RunY(RunStart(ii)+6)],'LineWidth',3);
            %                     plot(RunX(RunStart(ii):RunEnd(ii)),RunY(RunStart(ii):RunEnd(ii)),'k');
            %                     %scatter(RunX(RunStart(ii)),RunY(RunStart(ii)) ,3.5, HeadingStart1(ii),'filled');
            %
            %                     %colorbar
            %                 end
            %             end
            
            
            
 %--- get all bearing and heading angles for each left over run relative to
%%  preferred 02 isocline:------
            if ~isnan(RunX1{1})
                
                for R=1:length(RunX1)
                    
                    x=RunX1{R}(1:end);
                    y=RunY1{R}(1:end);
                    s=Tracks(1,T).Speed(RunStart(R):RunEnd(R));
                    s=s(1:3:end);
                    L=length(x);
                    bearingangles=NaN(1,1);
                    headingangles=NaN(1,1);
                    beelineX=NaN;
                    
                    
                    if F<1
                        CO2=2500;
                    else
                        CO2=10;
                    end
                    
                    for ii=1:L-1
                        beelineX=(CO2-(x(ii)));
                        %                     if F==5
                        %                          beelineX(ii)=(roi(1,1)-(y(ii)));
                        %                     end
                        beelineY=0;%(roi(2,2)+100-(x(ii)));  % alternatively bearing towards gas inlet
                        headingX=(x(ii+1))-(x(ii));
                        headingY=(y(ii+1))-(y(ii));
                        nenner=((beelineX*headingX)+(beelineY*headingY));
                        Brecher=(sqrt((beelineX).^2+(beelineY).^2))*(sqrt((headingX).^2+(headingY).^2));
                        cosinus_angle=nenner/Brecher;
                        bearingangles(ii)=(acos(cosinus_angle))*(360/(2*pi));
                        headingangles(ii) = atan2(headingY,headingX);
                    end
                    
                    cc=cc+1;
                    
                    Ctr_Bearing{T,R}=[bearingangles NaN ];
                    Ctr_dBearing{T,R}=[diff(bearingangles) NaN NaN ];
                    Ctr_Speed{T,R}=s;
                    
                end % end Run loop
                
            end
            
            %% --plot tracks--
%             if T>150 && T<190
%                 %
%                 hold on
%                 scatter(Tracks(1,T).SmoothX(1),Tracks(1,T).SmoothY(1),'dk')
%                 
%                 for R=1:length(RunX1)
%                     try
%                         scatter (RunX1{R}(1:end),RunY1{R}(1:end),30,Ctr_dBearing{T,R},'filled')
% %                         scatter(HeadRunsX{T,R},HeadRunsY{T,R},13,head_speed{T,R},'filled')
%                     catch;
%                     end
%                 end
%                 plot( Tracks(1,T).SmoothX, Tracks(1,T).SmoothY,'Color',[0.5 0.5 0.5]);
%                 title([T,F])
%                 colorbar
% 
%                 
%             end
            
        end
        
    end %end Tracks loop
    %%
    
    bearing{F}=Ctr_Bearing;
    dBearing{F}=Ctr_dBearing;
    speed_ctr{F}=Ctr_Speed;
     RunFrames{F}=RunFrames1;
    
end
%%

home=cd;
nd=(cd);
d= strfind(cd, '\');
name=nd(d(end)+1:end)

display('...save')

save(['runinfo_new' fname(1:end-42)  name] , 'bearing' ,'dBearing' ,'speed_ctr');



