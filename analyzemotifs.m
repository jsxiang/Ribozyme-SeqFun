addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/herrorbar
addpath ~/Documents/MATLAB/FACSseq/
addpath /Users/jx/Documents/MATLAB/FACSseq/FACSseqlib/
addpath /Users/jx/Documents/MFS/repo/
addpath ./cbrewer/
addpath ./subtightplot/

addpath ~/Documents/MATLAB/tSNE_matlab/
warning('off','all')



%%
r=load('RFS.mat');
rfs=r.rfs;


%% pairwise 
i=1;
for k=1:length(rfs)
percentile=5;
gaussian=0;
countthresh=10;
findPairwiseMu(rfs(k).alignedSdRD,rfs(k).loopII.mat,percentile,gaussian,{rfs(k).name_descriptive},countthresh)
set(gca,'linewidth',3)
set(gca,'fontsize',18)
c=colorbar;
ylabel(c,sprintf('Standardized log_1_0(RNA/DNA) 5th percentile'))
set(c,'fontsize',18)
set(c,'linewidth',1.5)

title('')
xlabel('position i')
ylabel('position j')
title(rfs(k).name_descriptive)
end

%% mutual information
% for k=1:length(rfs)
% percentile=5;
% gaussian=0;
% countthresh=10;
% findPairwiseMu(rfs(k).alignedSdRD,rfs(k).loopII.mat,percentile,gaussian,...
%{rfs(k).name_descriptive},countthresh)
% set(gca,'linewidth',3)
% set(gca,'fontsize',18)
% c=colorbar;
% ylabel(c,sprintf('Standardized log_1_0(RNA/DNA) 5th percentile'))
% set(c,'fontsize',18)
% set(c,'linewidth',1.5)
% 
% title('')
% xlabel('position i')
% ylabel('position j')
% title(rfs(k).name_descriptive)
% end

%%
k=5
[Y,I]=sort(rfs(k).alignedSdRD);
sortedseqmat=rfs(k).loopII.mat(I,:);
plotentropy=false;

% repeat the above but this time, instead of ranking the sequences, ...
% take fixed activity bins
numint=25;
minY=min(Y);
maxY=prctile(Y,99.9);
maxYi=length(Y)-1;
Yrange=Y(end-1)-minY;
setfig('histogram');clf


h=histogram(Y(1:maxYi),-4:0.25:2);

hedges=h.BinEdges;
hold on

yL=ylim;
yval=linspace(yL(1),yL(2));
plotedges=[22:-1:9 ];

for i=plotedges(2:end)
    xval=ones(1,100)*hedges(i);
    plot(xval,yval,'k:','linewidth',2.5)
end
ylabel('frequency')
title(sprintf('standardized log_1_0(RNA/DNA)'))
set(gca,'fontsize',22)
set(gca,'linewidth',1.5)

setfig('MI subplots hist2');clf
j=1;


for i=[plotedges]
    subtightplot(length(plotedges),1,j,0,0.05,0.3)
%     subtightplot(1,length(plotedges),j,0,0.3,0.05)
    if i==plotedges(end)
        xi=find(Y<hedges(i));
    elseif i==plotedges(1)
        xi=find(Y>hedges(i-1));
    else
        xi=find(Y<hedges(i)&Y>hedges(i-1));
    end
    findMutualInformation(sortedseqmat(xi,:),{'dummy'},false);
    L=ylim;
%     NumTicks=13;
%     set(gca,'YTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'Yticklabel',{'',1,'', 2,'', 3, '',4,'', 5, '',6,''})
%     
    
    NumTicks=11;
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    set(gca,'Yticklabel',{'',1,'', 2,'', 3, '',4,'', 5, ''})

    if i==16
        ylabel('position j')
    end
    
        xlabel('position i')
    set(gca,'fontsize',18)
    L=xlim;
    NumTicks=11;
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'xticklabel',{'',1,'', 2,'', 3, '',4,'', 5, ''})
    
%     NumTicks=13;
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'xticklabel',{'',1,'', 2,'', 3, '',4,'', 5, '',6,''})
%     
    
j=j+1;
end
% c=colorbar;
% ylabel(c,'mutual information')
% set(c,'fontsize',18)
% set(c,'linewidth',1.5)




%% output binned sequences for making seqlogo in R
for k=1:length(rfs)
    [Y,I]=sort(rfs(k).alignedSdRD);

sortedseq=rfs(k).loopII.seqs(I);
j=1;
for i=[plotedges]
    subtightplot(length(plotedges),1,j,0,0.05,0.3)
    if i==plotedges(end)
        xi=find(Y<hedges(i));
    elseif i==plotedges(1)
        xi=find(Y>hedges(i-1));
    else
        xi=find(Y<hedges(i)&Y>hedges(i-1));
    end
    if j<10
    filename=sprintf('%ssubseqs_0%d.txt',rfs(k).name,j);
    else
    filename=sprintf('%ssubseqs_%d.txt',rfs(k).name,j);
    end
    
    fid=fopen(filename,'w');
    for p=1:length(xi)
        if isempty(regexp(sortedseq{xi(p)},'w'))
        s=sortedseq{xi(p)};
        s=strrep(s,'a','A');
        s=strrep(s,'u','U');
        s=strrep(s,'c','C');
        s=strrep(s,'g','G');
        fprintf(fid,'%s\n',s);
        else
            continue;
        end
    end
    fclose(fid);
j=j+1;
end
end


%% compare N5 and N6
k=2;
loopIIlen=[];
for i=1:length(rfs(k).loopII.seqs)
    l2bases=regexprep(rfs(k).loopII.seqs{i},'-','');
    loopIIlen(i)=length(l2bases);
end
    
setfig('N5N6');clf
hold on
% histogram(rfs(k).alignedRDfold(loopIIlen==6),'binwidth',0.2)
% histogram(rfs(k).alignedRDfold(loopIIlen==5),'binwidth',0.2)
histogram(rfs(k).alignedSdRD(loopIIlen==6),'binwidth',0.2)
histogram(rfs(k).alignedSdRD(loopIIlen==5),'binwidth',0.2)

xlabel('Standardized log_1_0(RNA/DNA)(-ligand)')
ylabel('Frequency')
set(gca,'fontsize',24)
set(gca,'linewidth',2)
box on
title(rfs(k).name_descriptive)
legend('N6','N5','location','northeast')
xlim([-6 4])
setfig('N5N6 zoom');clf
hold on
% histogram(rfs(k).alignedRDfold(loopIIlen==6),'binwidth',0.2)
% histogram(rfs(k).alignedRDfold(loopIIlen==5),'binwidth',0.2)
histogram(rfs(k).alignedSdRD(loopIIlen==6),'binwidth',0.2)
histogram(rfs(k).alignedSdRD(loopIIlen==5),'binwidth',0.2)

% axis([-4.1 -2.7 -5 25])
set(gca,'fontsize',32)
set(gca,'linewidth',3)
box on

setfig('N5N6 fold');clf
hold on
histogram(10.^rfs(k).alignedRDfold(loopIIlen==6),'binwidth',0.1)
histogram(10.^rfs(k).alignedRDfold(loopIIlen==5),'binwidth',0.1)

xlabel('Fold (RNA/DNA)')
ylabel('Frequency')
set(gca,'fontsize',24)
set(gca,'linewidth',2)
box on
title(rfs(k).name_descriptive)
legend('N6','N5','location','northwest')

setfig('N5N6 fold zoom');clf
hold on
% histogram(rfs(k).alignedRDfold(loopIIlen==6),'binwidth',0.2)
% histogram(rfs(k).alignedRDfold(loopIIlen==5),'binwidth',0.2)
histogram(10.^rfs(k).alignedRDfold(loopIIlen==6),'binwidth',0.1)
histogram(10.^rfs(k).alignedRDfold(loopIIlen==5),'binwidth',0.1)


axis([1.8 2.8 -2 50])
set(gca,'fontsize',32)
set(gca,'linewidth',3)
box on








%% END







































