function [mI,pJoint,pointMI]=findMutualInformation(seqs,name,plotent)
% function calculates mutual information 
% input is an aligned array of sequences in 1,2,3,4
% output is a 2D vector of mutual information 
% MI=sum(Px,y(X,Y)*log2(Px,y(X,Y)/Px(X)Py(Y)))
% and joint probability Px,y(X,Y)
% load('MyColormaps','mycmap')
mycmap=cbrewer('seq','Blues',128);
% mycmap=mycmap(end:-1:1,:);
mycmap=[[1.00 1.00 1.00];mycmap];
% mycmap=mycmap(end:-1:1,:);
DNA={'A','U','C','G'};
l1l2seqs=seqs;

pJoint=zeros(length(DNA)*length(l1l2seqs(1,:)));
for i=1:length(DNA)
    for j=1:length(DNA)
        for k=1:length(l1l2seqs(1,:))
            for h=1:length(l1l2seqs(1,:))
                p=sum(l1l2seqs(:,k)==i&l1l2seqs(:,h)==j)/length(l1l2seqs(:,1));
                pJoint(4*(k-1)+i,4*(h-1)+j)=p;
                
            end
        end
    end
end
                
% try
%     figtitle=strcat(name,': joint probability');
% catch
%     figtitle='joint probability';
% end
% setfig(figtitle{1});clf
% imagesc(pJoint)
% c=colorbar;
% ylabel(c,'P(X,Y)')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% 
% mycmap=[[1.0 1.0 1.0];mycmap];
% colormap(mycmap)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('joint probability')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% 
% NumTicks = length(a)+1;
% 
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% 
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
% for i=1:length(mymajorgrids)
%     line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%     line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
% end
%     
    
mI=zeros(length(l1l2seqs(1,:)));
pointMI=zeros(size(pJoint));
try
[ent,pAll]=findEnt(seqs,name,plotent);
catch
    [ent,pAll]=findEnt(seqs,[],plotent);
end
for k=1:length(l1l2seqs(1,:))
    for h=1:length(l1l2seqs(1,:))
        all16mi=zeros(length(DNA));
        for i=1:length(DNA)
            for j=1:length(DNA)
                all16mi(i,j)=pJoint(4*(k-1)+i,4*(h-1)+j)*log2(pJoint(4*(k-1)+i,4*(h-1)+j)/pAll(i,k)/pAll(j,h));
                pointMI(4*(k-1)+i,4*(h-1)+j)=log2(pJoint(4*(k-1)+i,4*(h-1)+j)/pAll(i,k)/pAll(j,h));
            end
        end
        mI(k,h)=sum(sum(all16mi));
    end
end

% try
%     figtitle=strcat(name,': mutual information');
% catch
%     figtitle='mutual information';
% end
% setfig(figtitle{1});clf

% lc=length(mycmap(:,1));
% maxMI=0.25;
% if max(max(mI))>maxMI
%     fprintf('new max MI = %0.3f\n',max(max(mi)));
%     maxMI=max(max(mI));
% end
% max(max(mI));
% sc=round(max(max(mI))/maxMI*lc);
% mynewcmap=mycmap(1:sc,:);

imagesc(mI)
% colormap(mynewcmap)

colormap(mycmap)
set(gca,'fontsize',18)
set(gca,'linewidth',1.5)

% c=colorbar;
% ylabel(c,'mutual information')
% set(c,'fontsize',18)
% set(c,'linewidth',1.5)
% xlabel('position i')
% ylabel('position j')
% title(figtitle{1},'interpreter','none')


% try
%     figtitle=strcat(name,': all information');
% catch
%     figtitle='all information';
% end
% setfig(figtitle{1});clf
% imagesc(pointMI)
% c=colorbar;
% ylabel(c,'position specific information')
% set(c,'fontsize',12)
% set(c,'linewidth',1.5)
% 
% colormap(mycmap)
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% title('position specific information')
% L = get(gca,'XLim');
% xlabelnames=DNA;
% a=repmat({'','A','','T','','C','','G'},1,length(l1l2seqs(1,:)));
% 
% NumTicks = length(a)+1;
% 
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% set(gca,'xticklabel',{a{:},''})
% set(gca,'YTick',linspace(L(1),L(2),NumTicks))
% set(gca,'Yticklabel',{a{:},''})
% 
% mymajorgrids=0.5:4:(length(l1l2seqs(1,:))*4+0.5);
% 
% for i=1:length(mymajorgrids)
%     line([mymajorgrids(i) mymajorgrids(i)],ylim,'color','k')
%     line(xlim,[mymajorgrids(i) mymajorgrids(i)],'color','k')
% end

end