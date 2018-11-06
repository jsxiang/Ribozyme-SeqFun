addpath ~/Documents/robot/Matlab-Utilities/
addpath ~/Documents/MATLAB/
addpath ~/Documents/MATLAB/FACS/
addpath ~/Documents/MATLAB/herrorbar
addpath ~/Documents/MATLAB/FACSseq/
addpath /Users/jx/Documents/MATLAB/FACSseq/FACSseqlib/
addpath /Users/jx/Documents/MFS/repo/
addpath ~/Documents/MATLAB/cbrewer/cbrewer/cbrewer/
addpath ~/Documents/MATLAB/tSNE_matlab/
warning('off','all')

% rs=load('rslib.mat'); % load compiled RNAseq data
% rs2=load('rslib2.mat'); % load compiled RNAseq data
rs2=load('rslib2_new.mat'); % load compiled RNAseq data
rs2=rs2.rs;

%% analyze xanthine loop II motifs
k=2;
 
% top15xan=rs.lib(k).loopII.RD>prctile(rs.lib(k).loopII.RD,80) & rs.lib(k).loopII.RD<prctile(rs.lib(k).loopII.RD,90);
% bottom85xan=rs.lib(k).loopII.RD>prctile(rs.lib(k).loopII.RD,90) & rs.lib(k).loopII.RD<prctile(rs.lib(k).loopII.RD,100);
% [mi15,pJ15]=findMutualInformation(rs.lib(k).loopII.mat(top15xan,:),{'xan_<10'});
% [mi85,pJ85]=findMutualInformation(rs.lib(k).loopII.mat(bottom85xan,:),{'xan_90to100'});
% findDeltaMI(mi15,mi85,pJ15,pJ85,{'xanc

rs2.lib(k).loopII.RD=rs2.lib(k).alignedSdRD;
findPairwiseMu(rs2.lib(k).loopII.RD,rs2.lib(k).loopII.mat,5,0,{'xan'},10)
set(gca,'linewidth',3)
set(gca,'fontsize',16)
c=colorbar;
ylabel(c,sprintf('5th percentile \n Standardized log_1_0(RNA/DNA)'))
set(c,'fontsize',16)
set(c,'linewidth',1.5)

title('')
xlabel('position i')
ylabel('position j')

%%
figure(12)
c=colorbar;
ylabel(c,'Standardized log_1_0(RNA/DNA)')
set(c,'fontsize',20)

%%
addpath ~/Documents/MATLAB/subtightplot/subtightplot/
[Y,I]=sort(rs2.lib(k).loopII.RD);
sortedseqmat=rs2.lib(k).loopII.mat(I,:);
plotentropy=false;

numint=20;
xint=1:round(length(Y)/numint):length(Y);
setfig('MI subplots');clf
for i=2:length(xint)
%     subplot(1,length(xint)-1,i-1)
    subtightplot(1,length(xint)-1,i-1,0,0.4,0.1)
    [a,b,c]=findMutualInformation(sortedseqmat(xint(i-1):xint(i),:),{'dummy'},plotentropy);
    if i~=2
    set(gca,'Yticklabel',[])
    end
    if i==2
        ylabel('position j')
    end
    
    if i==round(numint/2)
        xlabel('position i')
    end
    set(gca,'fontsize',24)

end

c=colorbar;
ylabel(c,'mutual information')
set(c,'fontsize',18)
set(c,'linewidth',1.5)
%%
setfig('xan ranked');clf

area(1:length(Y),Y,'facecolor',[0.5,0.5,0.9],'edgecolor',[0.2,0.2,0.6],'linewidth',4)
xlim([1,length(Y)])
yL=ylim;
yval=linspace(yL(1),yL(2));

hold on
xint=1:round(length(Y)/(numint)):length(Y);

for i=xint
    xval=ones(1,100)*i;
    plot(xval,yval,'k:','linewidth',2)
    
end
xlabel('sorted xanthine library sequences')
ylabel(sprintf('standardized \nlog_1_0(RNA/DNA)'))

set(gca,'fontsize',20)
set(gca,'linewidth',2)
%% repeat the above but this time, instead of ranking the sequences, take fixed activity bins
numint=23;
minY=Y(1);
maxY=Y(end);
Yrange=maxY-minY;
setfig('histogram');clf

h=histogram(Y,'binwidth',Yrange/numint);

hedges=h.BinEdges;
hold on

yL=ylim;
yval=linspace(yL(1),yL(2));

for i=9:(length(hedges)-4)
    xval=ones(1,100)*hedges(i);
    plot(xval,yval,'k:','linewidth',2.5)
end
ylabel('frequency')
title(sprintf('standardized log_1_0(RNA/DNA)'))
set(gca,'fontsize',22)
set(gca,'linewidth',1.5)

setfig('MI subplots hist2');clf
j=1;
for i=(length(hedges)-3):-1:9
%     subplot(1,length(xint)-1,i-1)
%     subtightplot((length(hedges)-3)-8,1,j,0,0.05,0.3)
    subtightplot(1,(length(hedges)-3)-8,j,0,0.3,0.05)
    if i==9
        xi=find(Y<hedges(i));
    elseif i==(length(hedges)-3)
        xi=find(Y>hedges(i-1));
    else
        xi=find(Y<hedges(i)&Y>hedges(i-1));
    end
    findMutualInformation(sortedseqmat(xi,:),{'dummy'},false);
    L=ylim;
    NumTicks=13;
    set(gca,'YTick',linspace(L(1),L(2),NumTicks))
    set(gca,'Yticklabel',{'',1,'', 2,'', 3, '',4,'', 5, '',6,''})
%     if i~=9
%     set(gca,'Yticklabel',[])
%     end
    if i==14
        ylabel('position j')
    end
    
%     if j==(length(hedges)-3)
        xlabel('position i')
%     end
    set(gca,'fontsize',18)
    L=xlim;
    NumTicks=13;
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    set(gca,'xticklabel',{'',1,'', 2,'', 3, '',4,'', 5, '',6,''})
%     set(gca,'xticklabel',{1, 2, 3,4,5,6})
j=j+1;
end

c=colorbar;
ylabel(c,'mutual information')
set(c,'fontsize',18)
set(c,'linewidth',1.5)


%%
addpath ~/Documents/MATLAB/SeqLogoFig/
sortedseq=rs2.lib(k).loopII.seqs(I);
setfig('seqlogo subplots');clf
j=1;
for i=(length(hedges)-3):-1:9
%     subplot(1,length(xint)-1,i-1)
    subtightplot((length(hedges)-3)-8,1,j,0,0.05,0.3)
    if i==9
        xi=find(Y<hedges(i));
    elseif i==(length(hedges)-3)
        xi=find(Y>hedges(i-1));
    else
        xi=find(Y<hedges(i)&Y>hedges(i-1));
    end
    if j<10
    filename=sprintf('xansubseqs_0%d.txt',j);
    else
    filename=sprintf('xansubseqs_%d.txt',j);
    end
    
    fid=fopen(filename,'w');
    for p=1:length(xi)
        
        s=sortedseq{xi(p)};
        s=strrep(s,'a','A');
        s=strrep(s,'u','U');
        s=strrep(s,'c','C');
        s=strrep(s,'g','G');
        fprintf(fid,'%s\n',s);
    end
    fclose(fid);
%     findMutualInformation(sortedseqmat(xi,:),{'dummy'},false);
%     [npos,h]=SeqLogoFig(sortedseq(xi))
%     L=ylim;
%     NumTicks=13;
%     set(gca,'YTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'Yticklabel',{'',1,'', 2,'', 3, '',4,'', 5, '',6,''})
%     if i~=9
%     set(gca,'Yticklabel',[])
%     end
    if i==14
        ylabel('position j')
    end
    
%     if j==(length(hedges)-3)
        xlabel('position i')
%     end
    set(gca,'fontsize',18)
%     L=xlim;
%     NumTicks=13;
%     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
%     set(gca,'xticklabel',{'',1,'', 2,'', 3, '',4,'', 5, '',6,''})
%     set(gca,'xticklabel',{1, 2, 3,4,5,6})
j=j+1;
end
c=colorbar;
ylabel(c,'mutual information')
set(c,'fontsize',18)
set(c,'linewidth',1.5)











