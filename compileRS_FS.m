warning('off','all')

clear

tr=load('./data/theo_RSVI.mat');
tr=tr.theo_all;

tr20=load('./data/theo_all_min20.mat');
tr20=tr20.theo_all;

xr=load('./data/xan_RSVIII.mat');
xr=xr.xan_all_RSVIII;
xr=rmfield(xr,'fold');

cr=load('./data/cdGII_N5-6_RSVIIandRSIX.mat');
cr=cr.crad;

crI=load('./data/cdGI_RSIX.mat');
crI=crI.cdGI_all;

fr=load('./data/FA_pCS408_RSVIII.mat');
fr=fr.FA_all;

fr20=load('./data/FA20_RSVIII.mat');
fr20=fr20.pJX20FA_all;


tm=load('./data/theo_MFSVc.mat');
tm=tm.theo_mean_switch;

xm=load('./data/xan_MFSVe.mat');
xm=xm.xan_mean_switch;

cm=load('./data/cdG_MFSVc.mat');
cm=cm.cdG_mean_switch;

fm=load('./data/FA_MFSVb.mat');
fm=fm.FA_mean_switch;

%% compare absolute values of controls between libraries
tr.name='theo';
tr20.name='theo_mcnt20';
xr.name='xan';
cr.name='cdG';
crI.name='cdGI';
fr.name='FA';
fr20.name='FA20';


tr.name_descriptive='Theophylline';
tr20.name_descriptive='Theophylline (min count 20)';
xr.name_descriptive='Xanthine';
fr.name_descriptive='Folinic acid (GFP)';
fr20.name_descriptive='Folinic acid (mCherry)';

cr.name_descriptive='cyclic di-GMP II';
crI.name_descriptive='cyclic di-GMP I';
tr1=rmfield(tr,'allcounts');
tr2=rmfield(tr1,'allseqs');
tr20_1=rmfield(tr20,'allcounts');
tr20_2=rmfield(tr20_1,'allseqs');

rs(1)=tr2;
rs(2)=xr;
rs(3)=fr;

rs(4)=rmfield(cr,'fold');
rs(5)=crI;
rs(6)=tr20_2;
rs(7)=fr20;


for i=1:length(rs)
    fid=fopen(strcat(rs(i).name,'_seqs.fasta'),'w');
    for j=1:length(rs(i).seqs)
        if length(rs(i).seqs{j})>60; % filter out control ribozymes
        s=regexprep(rs(i).seqs{j},'T','U');
        fprintf(fid,'>seq%d\n',j);
        fprintf(fid,'%s\n',s);
        end
    end
    
    
    fclose(fid);
end

%% normalize and standardize data
% take a look at whether the data is normally distributed
setfig('lognormal');clf
for i=1:length(rs)
   subplot(length(rs),1,i)
=   xval=(rs(i).RDratio(:,1)-3);
   xvalstand=Standard(xval,1,mean(xval),std(xval));
   histogram(xvalstand,'binwidth',0.20)
xlim([-4,3])
   rs(i).standardRD=xvalstand;
end



%% perform local alignment linsi in shell
% for i=1:length(rs)
%     infilename=strcat(rs(i).name,'_seqs.fasta');
%     outfilename=strcat(rs(i).name,'_seqs_linsi.fasta');
%     
%     cmdstr=sprintf('/usr/local/bin/mafft-linsi %s > %s ',infilename,outfilename);
%     
%     fprintf('%s\n',cmdstr);
%     system(cmdstr);
%     
% end


%% parse aligned sequence files
for i=1:length(rs)
    
outfilename=strcat('./data/',rs(i).name,'_seqs_linsi.fasta');
fid=fopen(outfilename);
s=textscan(fid,'%s');
fclose(fid);

algnseq={};
seqidx=find(~cellfun('isempty',regexp(s{1},'>seq')));
for j=1:(length(seqidx)-1)
    aseq='';
    j1=seqidx(j);
    while (j1+1)<seqidx(j+1)
        aseq=strcat(aseq,s{1}{j1+1});
        j1=j1+1;
    end
    algnseq{end+1}=aseq;
    
end
aseq='';
for j=(j1+1):(length(s{1})-1)
    aseq=strcat(aseq,s{1}{j+1});
end
algnseq{end+1}=aseq;


algnseqFL={};
loopII.seqs={};
loopII.mat=[];

for j=1:length(algnseq)
    uppercasenodash=regexprep(algnseq{j},'u','T');
    uppercasenodash=regexprep(uppercasenodash,'t','T');
    uppercasenodash=regexprep(uppercasenodash,'a','A');
    uppercasenodash=regexprep(uppercasenodash,'g','G');
    uppercasenodash=regexprep(uppercasenodash,'c','C');
    uppercasenodash=regexprep(uppercasenodash,'-','');

    algnseqFL{end+1}=uppercasenodash;

    o=regexp(algnseq{j},'(cugaugaguc[a|c|g|u])([a|c|g|u|-]+)([a|c|g|u]gacgaaacagc)','tokens');
    try
    loopII.seqs{end+1}=o{1}{2};
    catch
        loopII.seqs{end+1}='naw';
    end    
    
end
loopII.lengths=cellfun('length',loopII.seqs);

for j=1:length(loopII.seqs)
if ~strcmp(loopII.seqs{j},'naw')
    loopIIvec=[];
    for k=1:length(loopII.seqs{j})
        loopIIvec(end+1)=regexp('aucg-',loopII.seqs{j}(k));
    end
else
    loopIIvec=zeros(1,max(loopII.lengths));
end

loopII.mat=[loopII.mat; loopIIvec];
end



[c,ia,ib]=intersect(algnseqFL,rs(i).seqs);

rs(i).alignedseq={};
rs(i).alignedseq=algnseq(ia);
rs(i).alignedSEQT=algnseqFL(ia);
rs(i).loopII.seqs=loopII.seqs(ia);
rs(i).loopII.mat=loopII.mat(ia,:);
rs(i).alignedSdRD=rs(i).standardRD(ib);
rs(i).alignedRDfold=(rs(i).RDratio(ib,2)-rs(i).RDratio(ib,1));
rs(i).loopII.lengths=loopII.lengths;

[c,ia,ib]=intersect(rs(i).alignedSEQT,rs(i).seqs);

rep1RD=[rs(i).minus.RDratio(rs(i).minus.combidx,1) rs(i).plus.RDratio(rs(i).plus.combidx,1)];
rep2RD=[rs(i).minus.RDratio(rs(i).minus.combidx,2) rs(i).plus.RDratio(rs(i).plus.combidx,2)];

SdRDrep1=Standard(rep1RD,1,mean(rep1RD(:)),std(rep1RD(:)));
SdRDrep2=Standard(rep2RD,1,mean(rep2RD(:)),std(rep2RD(:)));

    rs((i)).alignedRDrep1=zeros(length(rs((i)).alignedSEQT),2);
    rs((i)).alignedRDrep1(ia,:)=rep1RD(ib,:);
    rs((i)).alignedRDrep1(rs((i)).alignedRDrep1==0)=nan;
    
    rs((i)).alignedRDrep2=zeros(length(rs((i)).alignedSEQT),2);
    rs((i)).alignedRDrep2(ia,:)=rep2RD(ib,:);
    rs((i)).alignedRDrep2(rs((i)).alignedRDrep2==0)=nan;

rs(i).alignedSdRDrep1=zeros(length(rs(i).alignedSEQT),2);
rs(i).alignedSdRDrep1(ia,:)=SdRDrep1(ib,:);
rs(i).alignedSdRDrep1(rs(i).alignedSdRDrep1==0)=nan;

rs(i).alignedSdRDrep2=zeros(length(rs(i).alignedSEQT),2);
rs(i).alignedSdRDrep2(ia,:)=SdRDrep2(ib,:);
rs(i).alignedSdRDrep2(rs(i).alignedSdRDrep2==0)=nan;

rs(i).alignedSdRDfold1=zeros(length(rs(i).alignedSEQT),1);
rs(i).alignedSdRDfold1(ia)=[SdRDrep1(ib,2)-SdRDrep1(ib,1)];
rs(i).alignedSdRDfold1(rs(i).alignedSdRDfold1==0)=nan;

rs(i).alignedSdRDfold2=zeros(length(rs(i).alignedSEQT),1);
rs(i).alignedSdRDfold2(ia)=[SdRDrep2(ib,2)-SdRDrep2(ib,1)];
rs(i).alignedSdRDfold2(rs(i).alignedSdRDfold2==0)=nan;

rs(i).alignedSdRDidx=ib;

end

%% align FACS-seq data with RNA-seq

rfs=rs;

id=[1 2 7 4 6 ];

fs=[tm xm fm cm tm];

for i=1:length(id)
    [c,ia,ib]=intersect(rfs(id(i)).alignedSEQT,fs(i).seqs);
    
    murep1=[fs(i).minus.VYBmus1' fs(i).plus.VYBmus1'];
    murep2=[fs(i).minus.VYBmus2' fs(i).plus.VYBmus2'];
    
    SdMu1=Standard(murep1,1,mean(murep1(:)),std(murep1(:)));
    SdMu2=Standard(murep2,1,mean(murep2(:)),std(murep2(:)));
    
    rfs(id(i)).alignedmu1=zeros(length(rfs(id(i)).alignedSEQT),2);
    rfs(id(i)).alignedmu1(ia,:)=murep1(ib,:);
    rfs(id(i)).alignedmu1(rfs(id(i)).alignedmu1==0)=nan;
    
    rfs(id(i)).alignedmu2=zeros(length(rfs(id(i)).alignedSEQT),2);
    rfs(id(i)).alignedmu2(ia,:)=murep2(ib,:);
    rfs(id(i)).alignedmu2(rfs(id(i)).alignedmu2==0)=nan;
    
    rfs(id(i)).alignedmufold1=zeros(length(rfs(id(i)).alignedSEQT),1);
    rfs(id(i)).alignedmufold1(ia,:)=10.^(murep1(ib,2)-murep1(ib,1));
    rfs(id(i)).alignedmufold1(rfs(id(i)).alignedmufold1==0)=nan;
    
    rfs(id(i)).alignedmufold2=zeros(length(rfs(id(i)).alignedSEQT),1);
    rfs(id(i)).alignedmufold2(ia,:)=10.^(murep2(ib,2)-murep2(ib,1));
    rfs(id(i)).alignedmufold2(rfs(id(i)).alignedmufold2==0)=nan;
    
    rfs(id(i)).alignedSdMu1=zeros(length(rfs(id(i)).alignedSEQT),2);
    rfs(id(i)).alignedSdMu1(ia,:)=SdMu1(ib,:);
    rfs(id(i)).alignedSdMu1(rfs(id(i)).alignedSdMu1==0)=nan;
    
    rfs(id(i)).alignedSdMu2=zeros(length(rfs(id(i)).alignedSEQT),2);
    rfs(id(i)).alignedSdMu2(ia,:)=SdMu2(ib,:);
    rfs(id(i)).alignedSdMu2(rfs(id(i)).alignedSdMu2==0)=nan;
    
    rfs(id(i)).alignedSdMuFold1=zeros(length(rfs(id(i)).alignedSEQT),1);
    rfs(id(i)).alignedSdMuFold1(ia)=[SdMu1(ib,2)-SdMu1(ib,1)];
    rfs(id(i)).alignedSdMuFold1(rfs(id(i)).alignedSdMuFold1==0)=nan;
    
    rfs(id(i)).alignedSdMuFold2=zeros(length(rfs(id(i)).alignedSEQT),1);
    rfs(id(i)).alignedSdMuFold2(ia)=[SdMu2(ib,2)-SdMu2(ib,1)];
    rfs(id(i)).alignedSdMuFold2(rfs(id(i)).alignedSdMuFold2==0)=nan;
        
        
end


%% prep data into onehot
id=[1 2 7 4 6];


for i=1:length(rfs)
        positioncol={};
        rfs(i).onehot=[];
    for j=1:max(rfs(i).loopII.lengths)
        eachcol=[];
        for k=1:length(rfs(i).loopII.seqs)
            zvec=zeros(1,4);
            if rfs(i).loopII.mat(k,j)~=5 && rfs(i).loopII.mat(k,j)~=0
            zvec(rfs(i).loopII.mat(k,j))=1;
            end
            eachcol=[eachcol;zvec];

        end
        positioncol{j}=eachcol;
        rfs(i).onehot=[rfs(i).onehot eachcol];
    end

    rfs(i).alignedRDfold1=10.^(rfs(i).alignedRDrep1(:,2)-rfs(i).alignedRDrep1(:,1));
    rfs(i).alignedRDfold2=10.^(rfs(i).alignedRDrep2(:,2)-rfs(i).alignedRDrep2(:,1));
    rfs(i).onehotwlabel=[rfs(i).onehot rfs(i).alignedSdRDrep1 rfs(i).alignedSdRDrep2 rfs(i).alignedSdRDfold1 rfs(i).alignedSdRDfold2 ];
    rfs(i).RDmu=[rfs(i).alignedRDrep1 rfs(i).alignedRDrep2 rfs(i).alignedRDfold1 rfs(i).alignedRDfold2  ];
    if sum(i==id)>0
    rfs(i).onehotwlabel=[rfs(i).onehotwlabel rfs(i).alignedSdMu1 rfs(i).alignedSdMu2 rfs(i).alignedSdMuFold1 rfs(i).alignedSdMuFold2 ];
    
    rfs(i).RDmu=[rfs(i).RDmu rfs(i).alignedmu1 rfs(i).alignedmu2 rfs(i).alignedmufold1 rfs(i).alignedmufold2];
    end
end

%%
id=[1 2 7 4 6];
for i=1:length(rfs)
filename=sprintf('%s_onehotwALLlabels.csv',rs(i).name) ;
filename_nonSD=sprintf('%s_data.csv',rs(i).name) ;

header='';
formatcsv='';
for j=1:(max(rfs(i).loopII.lengths))
   header=strcat(header,sprintf('A%d,U%d,C%d,G%d,',j,j,j,j));
   formatcsv=strcat(formatcsv,'%d,%d,%d,%d,');
end

if sum(i==id)>0
    datalabel='RDrep1-,RDrep1+,RDrep2-,RDrep2+,RDfold1,RDfold2,mu-rep1,mu+rep1,mu-rep2,mu+rep2,mufold1,mufold2';
    header=strcat(header,datalabel);
    dataformatspec='%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n';
    formatcsv=strcat(formatcsv,dataformatspec);

    nonSDheader=sprintf('Sequence,%s',datalabel);
    nonSDformatcsv=strcat('%s,',dataformatspec);
else
% header=strcat(header,'RDrep1-,RDrep1+,RDrep2-,RDrep2+,RDfold1,RDfold2');
% formatcsv=strcat(formatcsv,'%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n');

    datalabel='RDrep1-,RDrep1+,RDrep2-,RDrep2+,RDfold1,RDfold2';
    header=strcat(header,datalabel);
    dataformatspec='%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n';
    formatcsv=strcat(formatcsv,dataformatspec);

    nonSDheader=sprintf('Sequence,%s',datalabel);
    nonSDformatcsv=strcat('%s,',dataformatspec);
    
end
fid=fopen(filename,'w');

fprintf(fid,'%s\n',header);
for j=1:length(rfs(i).onehotwlabel(:,1))
%    csvwrite(filename,rs.lib(i).onehotwlabel);
    fprintf(fid,formatcsv,rfs(i).onehotwlabel(j,:));

end
fclose(fid);

fid=fopen(filename_nonSD,'w');
fprintf(fid,'%s\n',nonSDheader);
for j=1:length(rfs(i).RDmu(:,1))
%    csvwrite(filename,rs.lib(i).onehotwlabel);
    fprintf(fid,nonSDformatcsv,rfs(i).alignedseq{j},rfs(i).RDmu(j,:));

end 
fclose(fid);

rfs(i).header=strsplit(header,',');

end

%% remove rows with nan for FA and cdG datasets
%%
id=[4 7];
for i=[4 7]
filename=sprintf('%s_onehotwmulabels.csv',rs(i).name) ;

header='';
formatcsv='';
for j=1:(max(rfs(i).loopII.lengths))
   header=strcat(header,sprintf('A%d,U%d,C%d,G%d,',j,j,j,j));
   formatcsv=strcat(formatcsv,'%d,%d,%d,%d,');
end

if sum(i==id)>0
header=strcat(header,'mu-rep1,mu+rep1,mu-rep2,mu+rep2,mufold1,mufold2');
formatcsv=strcat(formatcsv,'%1.6f,%1.6f,%1.6f,%1.6f,%1.6f,%1.6f\n');
end
fid=fopen(filename,'w');

fprintf(fid,'%s\n',header);
for j=1:length(rfs(i).onehotwlabel(:,1))
%    csvwrite(filename,rs.lib(i).onehotwlabel);
    if ~isnan(sum(rfs(i).onehotwlabel(j,31:36)))
    fprintf(fid,formatcsv,rfs(i).onehotwlabel(j,[1:24 31:36]));
    end

end
fclose(fid);


end

%%
% save('RFS.mat','rfs')
%%
k=6
RD1_pred=csvread(strcat('~/Documents/MFS/RS_MFS_all/runh2o/automl_360s/',rfs(k).name,'_RDrep1-_predicted.txt'),0,0);
RD2_pred=csvread(strcat('~/Documents/MFS/RS_MFS_all/runh2o/automl_360s/',rfs(k).name,'_RDrep2-_predicted.txt'),0,0);

setfig('h2o pred RD-');clf
hold on
i=find(~cellfun('isempty',regexp(rfs(k).header,'RDrep1-')));
j=find(~cellfun('isempty',regexp(rfs(k).header,'RDrep2-')));
plot(rfs(k).onehotwlabel(:,i),rfs(k).onehotwlabel(:,j),'.','color',[0.5 0.5 0.5],'markersize',15)

% plot(rfs(1).onehotwlabel(:,21),rfs(1).onehotwlabel(:,23),'.','color',[0.5 0.5 0.5],'markersize',15)
plot(RD1_pred,RD2_pred,'.','MarkerSize',15)

xlabel('RNA/DNA replicate 1')
ylabel('RNA/DNA replicate 2')
set(gca,'linewidth',2)
set(gca,'FontSize',24)
box on

mu1_pred=csvread(strcat('~/Documents/MFS/RS_MFS_all/runh2o/automl_360s/',rfs(k).name,'_mu-rep1_predicted.txt'),0,0);
mu2_pred=csvread(strcat('~/Documents/MFS/RS_MFS_all/runh2o/automl_360s/',rfs(k).name,'_mu-rep2_predicted.txt'),0,0);

setfig('h2o pred mu-');clf
hold on
m=find(~cellfun('isempty',regexp(rfs(k).header,'mu-rep1')));
n=find(~cellfun('isempty',regexp(rfs(k).header,'mu-rep2')));
plot(rfs(k).onehotwlabel(:,m),rfs(k).onehotwlabel(:,n),'.','color',[0.5 0.5 0.5],'markersize',15)

% plot(rfs(1).onehotwlabel(:,27),rfs(1).onehotwlabel(:,29),'.','color',[0.5 0.5 0.5],'markersize',15)
plot(mu1_pred,mu2_pred,'.','MarkerSize',15)

xlabel('mCherry/BFP replicate 1')
ylabel('mCherry/BFP replicate 2')
set(gca,'linewidth',2)
set(gca,'FontSize',24)
box on

setfig('h2o pred RD mu');clf
hold on


RD0=mean([rfs(k).onehotwlabel(:,i),rfs(k).onehotwlabel(:,j)],2,'omitnan');
mu0=mean([rfs(k).onehotwlabel(:,m),rfs(k).onehotwlabel(:,n)],2,'omitnan');

RD_pred=mean([RD1_pred,RD2_pred],2,'omitnan');
mu_pred=mean([mu1_pred,mu2_pred],2,'omitnan');

plot(RD0,mu0,'.','color',[0.5 0.5 0.5],'markersize',15)
plot(RD_pred(1:length(mu_pred)),mu_pred,'.','MarkerSize',15)
% plot1to1(RDfold_pred)
% [mdl0,rsq0]=plotlinfit(RD0,mu0)
% [mdl,rsq]=plotlinfit(RDfold_pred,mufold_pred)
xlabel('RNA/DNA')
ylabel('mCherry/BFP')
set(gca,'linewidth',2)
set(gca,'FontSize',24)
axis([-6 3 -4 2])
box on


%

%%





