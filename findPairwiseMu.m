function [pairwisemu,pairwisesigma] = findPairwiseMu(mus, seqs, percentile,gaussian,name,countthresh)
% function calculates pairwise mu from FACS-seq data
% input mus are previously determined
% input seqs are aligned sequences
% output gives the 2D vector of pairwise mu at 25th percentile
% if guassian is specified, a normfit is performed; otherwise it outputs the
% percentile-th value
% load('MyColormaps','mycmap')

mycmap=cbrewer('seq','Reds',28);
mycmap=mycmap(end:-1:1,:);
% mycmap=[[1.00 1.00 1.00];mycmap];

DNA={'A','U','C','G'};
l1l2mus=mus;
l1l2seqs=seqs;

pairwisemu=zeros(length(DNA)*length(l1l2seqs(1,:)));
pairwisecount=pairwisemu;

for i=1:length(DNA)
    for j=1:length(DNA)
        for k=1:length(l1l2seqs(1,:))
            for h=1:length(l1l2seqs(1,:))
                ind=l1l2seqs(:,k)==i&l1l2seqs(:,h)==j;
                m=l1l2mus(ind);
                s=sum(ind);
                if gaussian
                    [gmu,gsig]=normfit(m);
                    pairwisemu(4*(k-1)+i,4*(h-1)+j)=gmu;
                    pairwisesigma(4*(k-1)+i,4*(h-1)+j)=gsig;
                    pairwisecount(4*(k-1)+i,4*(h-1)+j)=s;
                else
                    pairwisemu(4*(k-1)+i,4*(h-1)+j)=prctile(m,percentile);
                    pairwisecount(4*(k-1)+i,4*(h-1)+j)=s;
                end
                
            end
        end
    end
end

countthreshmat=repmat(countthresh,length(l1l2seqs(1,:))*4,length(l1l2seqs(1,:))*4);

abovethresh=pairwisecount>=countthreshmat;
pairwisemu(isnan(pairwisemu))=(min(min(pairwisemu))-0.01);
pairwisemu=pairwisemu.*abovethresh;
try
    figtitle=strcat(name,': pairwise mu');
catch
    figtitle='pairwise mu';
end



setfig(figtitle{1});clf
imagesc(pairwisemu)
c=colorbar;
colormap(mycmap)
colobarlabel=sprintf('%0.0fth percentile \\mu',percentile);
ylabel(c,colobarlabel)
set(c,'fontsize',12)
set(c,'linewidth',2.5)

set(gca,'fontsize',10)
set(gca,'linewidth',1.5)
title(figtitle{1},'interpreter','none')
L = get(gca,'XLim');
xlabelnames=DNA;
a=repmat({'','A','','U','','C','','G'},1,length(l1l2seqs(1,:)));
colormap(mycmap)

NumTicks = length(a)+1;

set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'xticklabel',{a{:},''})
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'Yticklabel',{a{:},''})

mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
colormap(mycmap)

for i=1:length(mymajorgrids)
    line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k','linewidth',2)
    line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k','linewidth',2)
end

end