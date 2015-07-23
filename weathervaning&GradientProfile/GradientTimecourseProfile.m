figure('visible','on')
clear
home=cd;
geotypes=cell(8);
genotypes(1:8)={''};

cc=1;

% to edit:
tb=10; %averaging timebin
nP=4; % number of profiles

 
% average over e experiments for t timebins:
load HMnorm_batch
dim=NaN(1,3);
for i=1:length(HMnorm_batch)
dim(i,:)=size(HMnorm_batch{i});
end
dim=max(dim);
% clearvars -except directory home dd HMnorm_batch SpeedM_batch angSpeedM_batch cc genotypes
c=1;

HMoE={NaN};
t=1;
for c=1:tb:length(HMnorm_batch{1})-tb
    % it starts with a structure containing the normalized (percentage of
    % animals at a given timepoint in a givenbin of the arena)
    mHM=NaN(dim(1),dim(2));
    
    for e= 1:length(HMnorm_batch)
        mHM1=squeeze(nanmean(HMnorm_batch{1,e}(:,:,c:c+tb),3)); %averaging over tb 
        mHM(1:size(mHM1,1),1:size(mHM1,2),e)=mHM1;
    end

  %each cell: averaged data for e experiments and one timebin
  bi=find (mHM==Inf);
  mHM(bi)=NaN;
    HMoE{t}=mHM;
 
    t=t+1;
end
 nd=(cd);
 d= strfind(cd, '\');
name=nd(d(end)+1:end);


cd (home)

%% plot mean ocupancy for sectors over nP time bins
%now the data are organized by timepoint: each cell contains n arenas for
%one time point (values: %!!):
J=jet(length(HMoE));
cc=1;

for i=1:floor(length(HMoE)/nP):length(HMoE)
    disp(i)
sHM=nansum(HMoE{i},1); %!! summing along one dimension of arena(previous error:taking mean instead of the sum
sHM=squeeze(sHM);
meanOE=nanmean(sHM,2);% averaging over experiments
sorted=sort(meanOE);
semOE=nanstd(sHM,2)/sqrt(e);
hold on

 l=shadedErrorBar2([],meanOE,semOE,{'Color',J(i,:)},0);
 if ~isempty(l.patch)
 lh(cc)=l.patch;
 lg{cc}=num2str(i);
 set(gca,'XTick',[1:3:length(meanOE)])
 set(gca,'XTickLabel',(round([21:-(17/length(meanOE)*3.2):4]*10))/10);
 cc=cc+1;
 end
end

ylim([0 12])
% xlim([1 30])
title(name)
legend(lh,lg)
xlabel('% O2')
ylabel (' % worms')

