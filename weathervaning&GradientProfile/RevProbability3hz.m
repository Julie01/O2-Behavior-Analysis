
clear
warning off
% global rois
% cd('mat')
rate=3;
rate=round((rate*10/3)/10);
home=cd;
nd=(cd);
d= strfind(cd, '\');
name=nd(d(2):d(end));

files=dir('*als.mat');
nd=(cd);
d= strfind(cd, '\');
name=nd(d(2):d(end));


reversals={NaN};

%rois:
load roinew


% create gradient matrix
gradient=zeros(1060,2560);
gl=21:-17/(2500):4;
low=4;
high=21;

EC=1;

for i=1:size(gradient,1)
    gradient(i,1:30)=high;
    gradient(i,31:2531)=gl;
    gradient(i,2532:end)=low;
end

% cd ..\
for F=1:length(files)
    
    F
    clearvars -except F exp roi  home name files bearing dBearing speed_ctr...
        HeadSens HeadSpeed HeadSpeedNet reversals rate gradient sensPath
    cc=0;
    roiF=round(roi{F});
    cla
    fname=(files(F,1).name)
    load (files(F).name);
    %load([fname(1:end-4), '_anglesV8.mat']);
    Reversals={NaN};
    dCdT={NaN};
    
    pixToMM=0.0155;
    
    
    %% ---analyze Tracks----
    
    for T= 1:length(Tracks)
        % only tracks which start after establishment of the gradient
        % (ca frame 500)
        if  Tracks(1,T).Frames(end)>300 && Tracks(1,T).Frames(1)>100 
            
            if F<12
                TX=  (Tracks(1,T).SmoothY);
                TY=  (Tracks(1,T).SmoothX);
            else
                
                TX=  (Tracks(1,T).SmoothX);
                TY=  (Tracks(1,T).SmoothY);
            end
            
            
            %----exclude those which occur close to border:----
            bi1= find( TX<roiF(1,1)+200 | TX>roiF(2,1)-100);
            bi2= find( TY<roiF(2,2)+20 | TY>roiF(1,2)-20);
            bi= setdiff(bi1, bi2);
            bi=[bi,bi2];
            TX(bi)=NaN;
            TY(bi)=NaN;
            
            %kill runs which have unusual high angular speed
            %and don't cover 1 small worm travel distance (20 px):
            
            meanAngVelocity=(nanmean(abs(Tracks(1,T).AngSpeed)));
            if meanAngVelocity>30
                continue
            else
                d=NaN(1,1);
                ci=1;
                for iii=1:3:length(TY)
                    d(ci)=sqrt(((TX(1)-TX(iii)).^2)+((TY(1)-TY(iii)).^2));
                    ci=ci+1;
                end
                if max(d)<15
                    disp('no displacement')
                    continue
                end
            end
            
            %get O2 concentration at this position of path:
            try
                C_ind=sub2ind(size(gradient),TY(1:rate:end),TX(1:rate:end));
                nan_idx=isnan(C_ind);
                C_ind=C_ind(nan_idx~=1);
                
                %get O2 concentration at this position of path and reinsert Nans:
                sensoryPathNaN=(gradient(round(C_ind)));
                dCdT=NaN(1,length(nan_idx));
                dCdT(nan_idx~=1)=sensoryPathNaN;
                sensP{T}=dCdT;
            end
            
            % reversal vector:
            Revs=zeros(1,ceil(length(TX)/rate));
            Rev=(Tracks(1,T).polishedReversals);
            
            if ~isempty(Rev) & max(d)>=15;
                gi=find(Rev(:,4)>0);
                Rev=Rev(gi,1);
                Rev=round(Rev/rate);
                
                %Omegas=(Tracks(1,T).OmegaTrans(:,1));
                
                Revs(Rev)=1;
                
            end
            Reversals{T}=Revs;
            
            
            
            
            %--- get all bearing and heading angles for each left over run relative to
            %%  preferred 02 isocline:------
            if   max(d)>=15
                
                x=TX(1:rate:end);
                y=TY(1:rate:end);
                s=Tracks(1,T).Speed;
                s=s(1:rate:end);
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
                
                Ctr_Bearing{T}=[bearingangles NaN ];
                Ctr_dBearing{T}=[diff(bearingangles) NaN NaN ];
                Ctr_Speed{T}=s;
                
            end
            
        end
        
    end %end Tracks loop
    %%
    
    bearing{F}=Ctr_Bearing;
    dBearing{F}=Ctr_dBearing;
    speed_ctr{F}=Ctr_Speed;
    reversals{F}=Reversals;
    sensPath{F}=sensP;
    
end
%%

home=cd;
nd=(cd);
d= strfind(cd, '\');
name=nd(d(end)+1:end)

display('...save')

save(['revinfo_03hz_' fname(1:end-42)  name] , 'bearing' ,'dBearing' ,'speed_ctr','reversals','rate','sensPath');



