
files=dir('*als.mat');
load('roinew.mat');
roiF=roi;
roi={NaN(2,2)};


for F=1:length(files)
    
    load(files(F).name)
     roiFF=round(roiF{F});
    figure('Visible','on')
    
    for i= 1:2:length(Tracks)
        %plot the track and starting point:
        hold on       
        
        if length(Tracks(i).SmoothX)>150
            if F<12
                plot(Tracks(1,i).SmoothY,Tracks(1,i).SmoothX,'k')
            else
                plot(Tracks(1,i).SmoothX,Tracks(1,i).SmoothY,'k')
            end
            
        end       
        
    end
    title(F)
     axis image
     scatter(roiFF(1,1),roiFF(1,2))
scatter(roiFF(2,1),roiFF(2,2))
    pause (2)
   
%     roi{F}(1,:)=ginput(1);
%     roi{F}(2,:)=ginput(1);
   
     close
end

if ~isempty(files)
save roinew roi
end




