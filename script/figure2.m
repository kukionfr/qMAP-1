close all; clc; clear
runtag = 'figure2_250418';
%%
xlsname_datametrics='../data/skin_morphometric_data.xlsx';
xlsname_featureInfo='../data/feature_list.xlsx';
[num,~,all2]=xlsread(xlsname_datametrics,'by section all');
ccsi=strcmpi(all2(3,:),'SI'); %find columns with sample information
numlbl=all2(2:4,:); %column labels
% create data metrics
simat=num(:,ccsi); % sample info metrics; NaN for string 
simatc=all2(:,ccsi); %5th row and on contains data matrix
fmat=num(:,~ccsi); % feature metrics
silbl=numlbl(:,ccsi);
flbl=numlbl(:,~ccsi);
flbl=[num2cell(1:size(flbl,2));flbl]; % adding  identifyinig number in first row 

% sample info label; add one more column called sectionnum
silbl(3,end+1)={'SectionNum'};
silbl(2,end)={'SI'};
silbl(1,end)={nan}; 
silbl=[num2cell(1:size(silbl,2));silbl]; % adding identifyinig number
% ensemble patient data from sections.
cck=strcmpi(silbl(end,:),'patient'); %find patient ID column
up=unique(simat(:,cck));
up(isnan(up)) = [];
% up=up(isback&iscluerange);
plm=[]; pls=[]; pli=[]; plmed=[]; % patient level (PL) metrics, mean/std/info/median
plic={};
% iterate unique patient id
for kp=1:length(up)
    cc=simat(:,cck)==up(kp);
    pm=nanmean(fmat(cc,:),1);
    pmed=nanmedian(fmat(cc,:),1);
    ps=nanstd(fmat(cc,:),[],1);
    pi0=simat(cc,:);
    pi0=pi0(1,:);
    pi0c=simatc(cc,:);
    pi0c=pi0c(1,:);
    plm=[plm;pm]; % patient level mean (among sections)  
    plmed=[plmed;pmed]; % patient level median (among sections)
    pls=[pls;ps]; % patient level std (among secionts)
    pli=[pli;[pi0,sum(cc)]];
    plic=[plic;[pi0c,{sum(cc)}]];
end
%% get patietn sample info  
isMultiSec = cell2mat(plic(:,end))>1;

cc=1:size(flbl,2);
[cgender]=find(strcmpi(silbl(end,:),'gender'),1,'first');
ismale=strcmpi(plic(:,cgender),'male');
ccm=ismale; ccf=~ccm;

% select body part
[cbodypart]=find(strcmpi(silbl(end,:),'body part'),1,'first');
isback=strcmpi(plic(:,cbodypart),'back'); 

% select race
[crace]=find(strcmpi(silbl(end,:),'race'),1,'first');
iswhite=strcmpi(plic(:,crace),'white');

% select malignancy
[csamplingtype]=find(strcmpi(silbl(end,:),'Sample Behavior Type'),1,'first');
ismalignant=strcmpi(plic(:,csamplingtype),'Malignant'); 

% select biopsy type  
[csamplingproc]=find(strcmpi(silbl(end,:),'Sample Procedure Type'),1,'first');
isexcision=strcmpi(plic(:,csamplingproc),'Excis ion');

[cage]=find(strcmpi(silbl(end,:),'age')==1); %column id of age
age = pli(:,cage);

%% load tissue feature info 
[~,~,allf]=xlsread(xlsname_featureInfo);
allflbl=allf(1,:);
allf(1,:)=[];

cnonsense=strcmpi(allflbl,'nonsense');
isnonsense=cell2num(allf(:,cnonsense));

ctissuefeature=strcmpi(allflbl,'tissue_cell');
istissuefeature=strcmpi(allf(:,ctissuefeature),'tissue');
%% Correlation with aging and Cohen's D
cc=pli(:,end)>1;
plcv=abs(pls(:,:)./plm(:,:)); 
plcv(~cc,:)=nan;
% selecting samples used to calculate age correlation

ccsample = iswhite & isback; %& isMultiSec;
% ccsample = ~iswhite & isback;
% ccsample = ccsample & isold; %42 patients 

% setup criteron matrix and generate correlation matrix for each criteria
crimat{1}=plm>-inf; 
crimat{2}=plcv<0.5;  % CV below 0.5
crimat{3}=plcv<1;  % CV below 1

criname{1}='all';
criname{2}='CVless50';
criname{3}='CVless100';
    
lbltmp0={'N-all','N-male','N-female','rho-all', 'rho-male', 'rho-female',...
    'Pval-all','Pval-male','Pval-female','cohend-all','cohend-male','cohend-female',...
    'cohend-all-mal','cohend-male-mal','cohend-female-mal'};

plmuse=plmed;

cc=1:size(plmuse,2);

restmp0=[];
restmp0lbl={};
restmp0lbl2={};
tic;
ageoldth=60;
% Iterate criterion matrix
for kss=1:length(crimat)
    cef=0; cefm=0; ceff=0;
    p=0; pm=0; pf=0;
    Nm=0; Nf=0; Nall=0;
    cohendall=0;cohendallmal=0;
    cohendm=0;cohendmmal=0;
    cohendf=0;cohendfmal=0;
    % Iterate each feature column
    for kf=1:length(cc)
          % exclude nan value in the colum 
           tst= ccsample;
%            tst= ccsample & isMultiSec;
           tst= tst & ~isnan(plmuse(:,cc(kf)));
           tst= tst & crimat{kss}(:,cc(kf)); % adding set -criterion
           % can add more criteron such as stainning batch.. or exlcude
           % outlier
           Nall(kf)=sum(tst);
           Nm(kf)=sum(tst& ccm);
           Nf(kf)=sum(tst& ccf);

            if sum(tst>0)
                ageuse = pli(tst,cage);
                plmusek = plmuse(tst,cc(kf));
               [cef(kf),p(kf)]=corr(ageuse,plmusek,'type','spearman');
                isold0 = ageuse>=ageoldth; isyoung0 = ageuse<30;ismal0 = ismalignant(tst);
                matA = plmusek(isold0,:);matB = plmusek(isyoung0,:);
                cohendall(kf) = cohend(matA,matB);
                matA = plmusek(ismal0&~isyoung0,:);matB = plmusek(~ismal0&~isyoung0,:);
                cohendallmal(kf) = cohend(matA,matB);
            else
                cef(kf)=0;
                p(kf)=0;
                cohendall(kf)=0;
                cohendallmal(kf)=0;
            end

            if sum(tst & ccm>0)
                ageuse = pli(tst&ccm,cage);
                plmusek = plmuse(tst & ccm,cc(kf));
               [cefm(kf),pm(kf)]=corr(ageuse,plmusek,'type','spearman');
                isold0 = ageuse>=ageoldth; isyoung0 = ageuse<30;ismal0 = ismalignant(tst&ccm);
                matA = plmusek(isold0,:);matB = plmusek(isyoung0,:);
                cohendm(kf) = cohend(matA,matB);
                matA = plmusek(ismal0&~isyoung0,:);matB = plmusek(~ismal0&~isyoung0,:);
                cohendmmal(kf) = cohend(matA,matB);
            else
                cefm(kf)=0;
                pm(kf)=0;
                cohendm(kf)=0;
                cohendmmal(kf)=0;
            end

           if sum(tst&ccf>0)
               ageuse = pli(tst&ccf,cage);
               plmusek = plmuse(tst & ccf,cc(kf));
               [ceff(kf),pf(kf)]=corr(ageuse,plmusek,'type','spearman');
                isold0 = ageuse>=ageoldth; isyoung0 = ageuse<30;ismal0 = ismalignant(tst&ccf);
                matA = plmusek(isold0,:);matB = plmusek(isyoung0,:);
                cohendf(kf) = cohend(matA,matB);
                matA = plmusek(ismal0&~isyoung0,:);matB = plmusek(~ismal0&~isyoung0,:);
                cohendfmal(kf) = cohend(matA,matB);
           else
               ceff(kf)=0;
               pf(kf)=0;
               cohendf(kf)=0;
               cohendfmal(kf)=0;
           end            
    end
    restmp0=[restmp0, [Nall' Nm' Nf' cef' cefm'  ceff' p' pm' pf' cohendall' cohendm' cohendf' cohendallmal' cohendmmal' cohendfmal']]; %ceff and pf are zero
    restmp0lbl=[restmp0lbl,lbltmp0];
    restmp0lbl2=[restmp0lbl2,repmat(criname(kss),1,length(lbltmp0))];
end
restmp=[flbl(:,cc)' num2cell(restmp0)];
reslbltmp=cell(1,size(restmp,2));
reslbltmp(1,size(flbl,1)+1:end)=restmp0lbl2;
reslbltmp(2,size(flbl,1)+1:end)=restmp0lbl; 
if 1
    xlsname_corrmet=['../output/Feature Age correlation metrics',runtag,'.xlsx'];
    xlswrite(xlsname_corrmet,restmp,'correlation','A3');    
    xlswrite(xlsname_corrmet,reslbltmp,'correlation','A1');
end
toc;
%% Figure 2B - Sample cohort
ages = age(ccsample);
isMale = ccm(ccsample);

edges = 0:10:100;
ageGroups = discretize(ages, edges);

valid = ~isnan(ageGroups);
ageGroups = ageGroups(valid);
genderIdx = double(~isMale(valid)) + 1;  % Male = 1, Female = 2
counts = accumarray([ageGroups, genderIdx], 1, [numel(edges)-1, 2]);

bar(counts, 'grouped');
xticklabels(arrayfun(@(a,b) sprintf('%d-%d', a+1, b), edges(1:end-1), edges(2:end), 'UniformOutput', false));
legend('Male', 'Female');bjff3;
saveas(gcf, '../output/figure2b.png');
saveas(gcf, '../output/figure2b.fig');
close all;
%% Figure 2C - Volcano plot for feature selection
fall=cell2mat(restmp(:,[8 11 14])); %8:rho 11:p-val 14:cohend

rho = fall(:,1);
pval = fall(:,2);
cd = fall(:,3);

cc = abs(cd)>0.8 & abs(rho)>0.3 & pval<0.05; 
figure;plot(rho(~cc),abs(cd(~cc)),'k.');hold on;plot(rho(cc),abs(cd(cc)),'r.');bjff3;
xticks([-1:0.5:1]);xlim([-1 1]);ylim([0 2.5]);grid on;
saveas(gcf, '../output/figure2c.png');
saveas(gcf, '../output/figure2c.fig');
close all;
%% Figure 2D - known aging feature
ccf = abs(fall(:,1))>0.3 & fall(:,2)<0.05 & abs(fall(:,3))>0.8 & ~isnonsense;
ageuse = age(ccsample);
plmuse = plmed(ccsample,:);
figure;plot(ageuse,plmuse(:,116),'k+');bjff3;
saveas(gcf, '../output/figure2d_1.png');
saveas(gcf, '../output/figure2d_1.fig');
close all
figure;plot(ageuse,plmuse(:,76),'k+');bjff3;
saveas(gcf, '../output/figure2d_2.png');
saveas(gcf, '../output/figure2d_2.fig');
close all
%% Figure 2E - novel aging feature
figure;plot(ageuse,plmuse(:,36),'k+');bjff3;
saveas(gcf, '../output/figure2e_1.png');
saveas(gcf, '../output/figure2e_1.fig');
close all
figure;plot(ageuse,plmuse(:,170),'k+');bjff3;
saveas(gcf, '../output/figure2e_2.png');
saveas(gcf, '../output/figure2e_2.fig');
close all
%% Figure 2F - clustering

% preprocess data t
plmuse(isnan(plmuse))=0; %replace inf to zero 
[agesorted,I]=sort(ageuse); %sort with age just in case
morph_data = plmuse(I,ccf);
morph_data_norm =zscore(morph_data);
%pairwise correlation
[rho,pval] = corr(morph_data_norm,'Type','Spearman','Rows','pairwise');

cg2 = clustergram(rho,'Colormap',redbluecmap,'RowPdist', 'cityblock',...
                             'ColumnPdist', 'cityblock');
set(cg2,'Dendrogram',18)

% Note*** g1~8 need to be manually defined within clustergram 
% right click on tree node > export to workspace 
manual_select=0;
if manual_select 
    gs = {g1, g2, g3, g4, g5, g6, g7, g8};
    allTables = [];
    
    for i = 1:length(gs)
        rowLabels = gs{i}.RowLabels;
        featureID = str2double(strtrim(rowLabels)); 
        groupID = i * ones(length(featureID), 1);  % Assign group i
        T = table(featureID, groupID);
        allTables = [allTables; T];  % Append to the full table
    end
    allTables = sortrows(allTables, 'featureID');
    featID = find(ccf);
    T2 = cell2table(allf(featID,:));
    T2.Properties.VariableNames = allflbl;
    allTables = [allTables T2];
    allTables.cohenD=cd(ccf)';
    writetable(allTables, '../output/clustered_aging_features.csv');
end
%% output data for figure 3 and 4.
save('../output/aging_feature.mat', 'morph_data_norm','plmed','ccf','ccsample');
save('../output/sample_ages.mat', 'ageuse');
ismale = ccm(ccsample);
save('../output/sample_sex.mat','ismale');
flbls = flbl(:,ccf);
save('../output/feature_label.mat','flbl','flbls');
fstat = fall(ccf,:);
save('../output/feature_stats.mat','fstat');
