save('../output/aging_feature.mat', 'morph_data_norm','plmed','ccf','ccsample');
load('../output/sample_ages.mat', 'ageuse');
load('../output/sample_sex.mat','ismale');
load('../output/feature_label.mat','flbl','flbls');
load('../output/feature_stats.mat','fstat');
%% FIG 3A: plot rho vs univaraite MAE of correlated features
runtag = '250419';
morph_data = plmed(ccsample,ccf);
pmresuf = univariateglm(morph_data,ageuse,runtag,flbls,ismale);
sprintf('MAE mean %f std %f',nanmean(pmresuf(:,1)),nanmean(pmresuf(:,2)))

MAEint = pmresuf(:,1);
rho = fstat(:,1);
figure;
p = plot(rho,MAEint,'k.');bjff3;
p.MarkerSize = 10;
xlim([-0.9 0.9]);ylim([12 18])
saveas(gcf, '../output/figure3a.png');
saveas(gcf, '../output/figure3a.fig');
close all

[tmp,edges]=histcounts(MAEint,[12:.1:18]);
bar(tmp);bjff3;ylim([0 14])
set(gca, 'XDir','reverse');xticks([]);set(gca,'YAxisLocation','right')
set(gcf, 'Position',[100 100 300 100])
%% FIG 3B: 
[minVal, minIdx] = min(MAEint);
% disp(flbl(:,minIdx));
best_uni = morph_data(:,minIdx);
figure;plot(ageuse,best_uni,'k+');bjff3;
saveas(gcf, '../output/figure3b_1.png');
saveas(gcf, '../output/figure3b_1.fig');
close all

x0=best_uni;
y0=ageuse;
ypreds=[];yts=[];
parfor kcv=1:length(y0) %leave-one-out model for each feature
    p=[1:kcv-1,kcv+1:length(y0),kcv];
    trainnum=length(p)-1;
    trainid=p(1:trainnum);
    testid=p(trainnum+1:end);
    x=x0(trainid,:);
    y=y0(trainid);
    xt=x0(testid,:);
    yt=y0(testid);
    % glm
    mdl = fitglm(x,y);
    ypred = predict(mdl,xt);
    ypreds=[ypreds;ypred];
    yts =[yts,yt];
end 
figure;plot(yts,ypreds,'k+');bjff3;
xlim([0 100]);ylim([0 100]);
xticks([0:20:100]);yticks([0:20:100]);
saveas(gcf, '../output/figure3b_2.png');
saveas(gcf, '../output/figure3b_2.fig');
close all
%% FIG 3C
unierr = MAEint;
clustergroup = readtable('../output/clustered_aging_features.csv');
clgroup=clustergroup.groupID;
clgroup(isnan(clgroup))=2;
errors = cell(1,8);
for i = 1:8
    e = unierr(clgroup == i);
    errors{i} = e(~isnan(e));  % remove NaNs
end

violin_dots(errors, 1, [.8 .8 .8]);bjff3;
set(gcf,'position',[100 100 720 300]);
xticks([1:8]);xlim([0,9]);
ylim([10,20]);yticks([10:2:20]);grid on;
saveas(gcf, '../output/figure3c.png');
saveas(gcf, '../output/figure3c.fig');
close all
%% FIG 3D-F setup
% Setup
ogLUT = find(ccf);
runtag = '250419';
vnew = 1:sum(ccf);
C = nchoosek(vnew, 2);  % all bivariate feature combinations
nPairs = size(C, 1);
nSamples = length(y0);

% Preallocate results
pmresuf_all = zeros(nPairs, 7);  % [varID1 varID2 group1 group2 MAE permMAE stdErr]

% ---------------------- PROGRESS BAR SETUP ----------------------
% Create waitbar and tracker
hWaitbar = waitbar(0, 'Running pairwise GLMs...');
tracker = ProgressTracker(nPairs, hWaitbar);
dq = parallel.pool.DataQueue;
afterEach(dq, @(~) tracker.increment());
% ----------------------------------------------------------------

% Run parallel loop over feature pairs
parfor kff = 1:nPairs
    corfeatids = C(kff, :);  % current feature pair
    x0 = morph_data_norm(:, corfeatids);  % get features
    yerrs = zeros(nSamples, 1);
    randyerrs = zeros(nSamples, 1);
    yts = zeros(nSamples, 1);
    ypreds = zeros(nSamples, 1);

    for kcv = 1:nSamples
        train_idx = setdiff(1:nSamples, kcv);
        test_idx = kcv;

        x_train = x0(train_idx, :);
        y_train = y0(train_idx);
        x_test = x0(test_idx, :);
        y_test = y0(test_idx);

        % LM model
        mdl = fitlm(x_train, y_train);
        y_pred = predict(mdl, x_test);
        yerrs(kcv) = abs(y_pred - y_test);

        % Permuted LM
        y_train_perm = y_train(randperm(length(y_train)));
        mdl_perm = fitlm(x_train, y_train_perm);
        y_pred_perm = predict(mdl_perm, x_test);
        randyerrs(kcv) = abs(y_pred_perm - y_test);

        % Optional prediction storage
        yts(kcv) = y_test;
        ypreds(kcv) = y_pred;
    end

    % Extract original IDs and cluster groups
    ogfeatids = ogLUT(corfeatids)';
    clgroups = clgroup(corfeatids)';
    if length(clgroups) < 2
        clgroups = [clgroups, 0];
    end

    % Save result
    pmresuf_all(kff, :) = [ogfeatids, clgroups, ...
        mean(yerrs), mean(randyerrs), std(yerrs)];

    % Update progress bar
    send(dq, []);
end

% Close the waitbar
close(hWaitbar);

% Sort results by MAE
pmresuf = sortrows(pmresuf_all, 5);

% Prepare table for Excel
ranks = num2cell((1:size(pmresuf, 1))');

labelA = flbl(2:end, pmresuf(:, 1))';
labelB = flbl(2:end, pmresuf(:, 2))';


% Preallocate output array
longArray=pmresuf(:, 1);
indexArray = zeros(size(longArray));
[~, idx] = ismember(longArray, ogLUT);
indexArray = idx;
unierrA = unierr(indexArray);

longArray=pmresuf(:, 2);
indexArray = zeros(size(longArray));
[~, idx] = ismember(longArray, ogLUT);
indexArray = idx;
unierrB = unierr(indexArray);
unierrM=mean([unierrA,unierrB],2);

bipred = [ranks, labelA, labelB, num2cell(pmresuf),num2cell(unierrM)];

% Add headers
headers = {'rank', 'tissue component A', 'feature group A', 'variable name A', ...
           'tissue component B', 'feature group B', 'variable name B', ...
           'variable id A', 'variable id B', 'clgroup A', 'clgroup B', ...
           'MAE', 'permMAE', 'stdev_err','mean_univariate_error'};
bipred = [headers; bipred];

% Save to Excel
xlsname_bipred = ['../output/Bivariate Predicting power', runtag, '.xlsx'];
writecell(bipred, xlsname_bipred);
%%  FIG 3D
figure;
scale='log';
gradenum=100;
gradelimit=100;              
DensitySCPlotW(pmresuf(:, 5),unierrM,gradenum,gradelimit,scale);bjff3;axis equal;
xlim([13 18]);ylim([13 18])
hold on; plot([13:18],[13:18],'k--')
xticks([13:2:18]);yticks([13:2:18]);

saveas(gcf, '../output/figure3d.png');
saveas(gcf, '../output/figure3d.fig');
close all
%% FIG 3E: Best bivariate pair
% Sample data
y1 = plmed(ccsample,pmresuf(1,1));   
y2 = plmed(ccsample,pmresuf(1,2));  
x = ageuse;           % Shared x-axis
% Plot all three arrays
figure;
yyaxis left
plot(x, y1, 'r+');
ylabel('spinosum thickness');

yyaxis right
plot(x, y2, 'bs');
ylabel('inflammatory cell density');
xticks([0:20:100]);bjff3;

saveas(gcf, '../output/figure3e_1.png');
saveas(gcf, '../output/figure3e_1.fig');
close all
%% FIG 3E-2

corfeatids = pmresuf(1,1:2);  % current feature pair
[~, idx] = ismember(corfeatids, ogLUT);

    x0 = morph_data_norm(:, idx);  % get features
    yerrs = zeros(nSamples, 1);
    randyerrs = zeros(nSamples, 1);
    yts = zeros(nSamples, 1);
    ypreds = zeros(nSamples, 1);

    for kcv = 1:nSamples
        train_idx = setdiff(1:nSamples, kcv);
        test_idx = kcv;

        x_train = x0(train_idx, :);
        y_train = y0(train_idx);
        x_test = x0(test_idx, :);
        y_test = y0(test_idx);

        % LM model
        mdl = fitlm(x_train, y_train);
        y_pred = predict(mdl, x_test);
        yerrs(kcv) = abs(y_pred - y_test);

        % Permuted LM
        y_train_perm = y_train(randperm(length(y_train)));
        mdl_perm = fitlm(x_train, y_train_perm);
        y_pred_perm = predict(mdl_perm, x_test);
        randyerrs(kcv) = abs(y_pred_perm - y_test);

        % Optional prediction storage
        yts(kcv) = y_test;
        ypreds(kcv) = y_pred;
    end
%Bivariate MAE
figure;plot(yts,ypreds,'k.');hold on;
plot(1:100,1:100,'k--');bjff3;
xticks([0:20:100]);xticklabels([]);
yticks([0:20:100]);yticklabels([]);
saveas(gcf, '../output/figure3e_2.png');
saveas(gcf, '../output/figure3e_2.fig');
close all
%% FIG 3F: plot bivariate pairs by clgroup
clpair = pmresuf(1:100,3:4);
%count occurence of feature group and core process pair
[C,ia,ic] = unique(clpair,'row');
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

value_counts(value_counts(:,2)==0,:)=[];
Data = zeros(max(value_counts(:,1)),max(value_counts(:,2)));
for i = 1:size(value_counts,1)
    Data(value_counts(i,1),value_counts(i,2))=value_counts(i,3);
end
%Data = 8 x8 matrix  
myColorMap = zeros(length(Data),3);
% myColorMap = lines(length(Data));
nodeNames = {'', '', '', '', '', '', '', ''};
circularGraph(Data,'Colormap',myColorMap,'Label',nodeNames);
saveas(gcf, '../output/figure3f.png');
saveas(gcf, '../output/figure3f.fig');
close all
%% FIG 3G: 
tmp=zscore(plmed(ccsample,ccf));
tmp=plmed(ccsample,ccf);
tmp(isnan(tmp))=0;
morph_data_norm=zscore(tmp); %Nfeat=109

[coeff,score,latent]=pca(morph_data_norm);
vcum=cumsum(latent)/sum(latent);
v95num=find(vcum-0.95>0,1,'first');
x0=score(:,1:v95num);   %Nfeat=30

imagesc(x0);

saveas(gcf, '../output/figure3g_1.png');
saveas(gcf, '../output/figure3g_1.fig');
close all
%% FIG 3G: PC number
figure;plot(1:98,vcum,'k-');hold on;xline(v95num, ':k');yline(0.95, ':k');plot(v95num,0.95,'ko');
bjff3;xticks([0,v95num,100]);ylim([0 1]);yticks([0,0.5,0.95,1]);box off;  
saveas(gcf, '../output/figure3g_2.png');
saveas(gcf, '../output/figure3g_2.fig');
close all
%% FIG3G: PC1 vs PC2 for age groups
ageyoung = ageuse<30;
agemiddle = 30<=ageuse & ageuse<=60;
ageold = 60<ageuse;
[sum(ageyoung),sum(agemiddle),sum(ageold)]
figure;plot(x0(ageyoung,1),x0(ageyoung,2),'r.');hold on
plot(x0(agemiddle,1),x0(agemiddle,2),'g.');
plot(x0(ageold,1),x0(ageold,2),'b.');bjff3;
xlim([-12 12]);ylim([-12 12]);yticks([-10 0 10]);
saveas(gcf, '../output/figure3g_3.png');
saveas(gcf, '../output/figure3g_3.fig');
close all
%% FIG3H
[maes,errmat,y_glm,y_las,y_svm,y_rtree,randcorr] = agingmodel(x0,ageuse);
figure;violin_dots({abs(errmat(:,2)),abs(errmat(:,4)),abs(errmat(:,6))},1,[.8 .8 .8]);
ylim([0 40]);xticklabels([]);bjff3;
saveas(gcf, '../output/figure3h_1.png');
saveas(gcf, '../output/figure3h_1.fig');
close all
%% FIG3H
figure;plot(ageuse,y_svm,'k.');xticks([0:20:100]);xticklabels([]);
yticks([0:20:100]);yticklabels([]);bjff3;
saveas(gcf, '../output/figure3h_2.png');
saveas(gcf, '../output/figure3h_2.fig');
close all