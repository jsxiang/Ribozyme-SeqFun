function [ent,pAll]=findEnt(seqs,name,plotent)
% the input seqs is an alignment of sequences
% in the case of FACS-seqs they are the loops of the same lengths
% output figure shows the fraction of each base at each position
% height of the bars correspond to the level of entropy at each position
pA=sum(seqs==1)/length(seqs(:,1));
pT=sum(seqs==2)/length(seqs(:,1));
pC=sum(seqs==3)/length(seqs(:,1));
pG=sum(seqs==4)/length(seqs(:,1));

pAll=[pA;pT;pC;pG];
ent=2+pA.*log2(pA)+pT.*log2(pT)+pC.*log2(pC)+pG.*log2(pG);

addpath ~/Documents/MATLAB/cbrewer/cbrewer/cbrewer

try
    figtitle=strcat(name,': entropy');
catch
    figtitle={'entropy'};
end

if plotent
    mycmap=cbrewer('qual','Pastel1',8);

setfig(figtitle{1});clf
colormap(mycmap)

bar([pAll.*repmat(ent,4,1)]','Stacked','FaceAlpha',0.75)

title(figtitle{1},'interpreter','none')


legend('A','U','C','G','Location','Best')
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
ylabel('entropy')
xlabel('position')
end
end