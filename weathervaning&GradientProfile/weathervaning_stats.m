
folders=dir('*2*');

allTB=NaN(34,4);
sf={'control','gradient'};
cc=1;
  

for GT= 1:2
      
    cd(folders(GT).name)
    
    for i=1:length(sf)
    
    cd(sf{i})
    
    load turningbias
    
    allTB(1:length(mT),cc)=mT(:,1);
    cc=cc+1;
    
    cd ..\
    end
    
    cd ..\
    
end

figure
hold on
bar (nanmean(allTB))
errorb(nanmean(allTB),nanstd(allTB)/sqrt(34))
ylim([-0.2 0.05])
ylabel ('turning bias')
title ('central bins 40-150')
set(gca,'XTickLabel',sf);
set(gca,'XTick', 1:4)

for ii= [1 3]
[h(ii), p(ii)]=nanranksum(allTB(:,ii),allTB(:,ii+1));
text(ii+0.3,0.04,num2str(h(ii)))
end

cc=1
for ii= [1 2]
[h(cc), p(cc)]=nanranksum(allTB(:,ii),allTB(:,ii+2));
text(cc,-0.18,num2str(h(cc)))
cc=cc+1;
end






    
    