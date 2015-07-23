figure('visible','on')
clear
home=cd;

% to edit:

nP=3; % number of profiles

files= dir('*HMnorm*');
HMnorm_batch=[];
HMoE=cell(1,30);
mHM=NaN(1,37);
E=1;
EC=0;
e=0;
for F= 1:3%length(files)
    
load (files(F).name);
disp(files(F).name)
%     set1=set1.HMnorm_batch;
%  HMnorm_batch=[HMnorm_batch,set1];   
 

% average 10 timepoints:
tb=round(length(HMnorm_batch{1})/30); %averaging timebin
dim=NaN(1,3);
for i=1:length(HMnorm_batch)
    dim(i,:)=size(HMnorm_batch{i});
end
dim=max(dim);
% clearvars -except directory home dd HMnorm_batch SpeedM_batch angSpeedM_batch cc genotypes
c=1;

t=1;

for c=1:tb:length(HMnorm_batch{1})-tb
    % it starts with a structure containing the normalized (percentage of
    % animals at a given timepoint in a givenbin of the arena)
   
    %each cell: averaged data for e experiments and one timebin
     mHM=HMoE{t};
    E=E-e;
    for e= 1:length(HMnorm_batch)
        mHM1=squeeze(nanmean(HMnorm_batch{1,e}(:,:,c:c+tb),3)); %averaging over tb
        %interpolate for equal point number in X dimension:
        Y=nansum(mHM1);
        Xup=interp1(1:length(Y),Y,1:length(Y)/37:length(Y));
        mHM(1,1:length(Xup),E)=Xup;
        E=E+1;
    end
    
    bi=find (mHM==Inf);
    mHM(bi)=NaN;
    HMoE{t}=mHM;
    t=t+1;
    
end
e=0;

end % end files loop
    
 nd=(cd);
 d= strfind(cd, '\');
name=nd(d(end)+1:end) ;
nameL=[nd(d(end-1)+1:d(end)-1) ': ' nd(d(end)+1:end)];


cd (home)

%% plot mean ocupancy for sectors over nP time bins
%now the data are organized by timepoint: each cell contains n arenas for
%one time point (values: %!!):
J=jet(length(HMoE)*1.3);
cc=1;
nP=nP-1;

for i=1:floor(length(HMoE)/nP)-1:length(HMoE)
    disp(i)
sHM=(HMoE{i}); %!! summing along one dimension of arena(previous error:taking mean instead of the sum
sHM=squeeze(sHM);
meanOE=nanmean(sHM,2);% averaging over experiments
semOE=nanstd(sHM,2)/sqrt(E);
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
title(nameL)
legend(lh,lg)
xlabel('% O2')
ylabel (' % worms')

saveas(gca, [name '_profile.fig'])

