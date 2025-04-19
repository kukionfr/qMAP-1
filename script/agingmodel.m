function [maes,errmat,y_glm,y_las,y_svm,y_rtree,randcorr] = agingmodel(x0,y0)
    tic;
    errmat=[];y_glm=[];y_las=[];y_svm=[];y_rtree=[];%predictker=[];
    parfor kcv=1:length(y0)
        p=[1:kcv-1,kcv+1:length(y0),kcv];
        trainnum=length(p)-1;
        trainid=p(1:trainnum);
        testid=p(trainnum+1:end);
        x=x0(trainid,:);
        y=y0(trainid);
        randy = y(randperm(length(y)));
        xt=x0(testid,:);
        yt=y0(testid);
        % glm
        mdl = fitglm(x,y);
        ypred = predict(mdl,xt);
        yerr=(ypred-yt);
        randmdl = fitglm(x,randy);
        randypred = predict(randmdl,xt);
        randyerr=(randypred-yt);
        % lasso
        lambda = 1e-03;
        [B, FitInfo] = lasso(x,y,'Alpha',0.5,'CV',10);
        idxLambda1SE = FitInfo.Index1SE;
        idxLambda1SE = FitInfo.IndexMinMSE;
        coef = B(:,idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
        ypred2 = xt*coef + coef0;
        yerr2= (ypred2-yt);
        [B, FitInfo] = lasso(x,randy,'Alpha',0.5,'CV',10);
        idxLambda1SE = FitInfo.Index1SE;
        idxLambda1SE = FitInfo.IndexMinMSE;
        coef = B(:,idxLambda1SE);
        coef0 = FitInfo.Intercept(idxLambda1SE);
        randypred2 = xt*coef + coef0;
        randyerr2= (randypred2-yt);
        % svm
        svmmdl = fitrsvm(x,y,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
        ypred3 = predict(svmmdl,xt);
        yerr3=(ypred3-yt);
        randsvmmdl = fitrsvm(x,randy,'KernelFunction','linear','Epsilon',0.1,'Alpha',ones(size(y,1),1)*0.001 );
        randypred3 = predict(randsvmmdl,xt);
        randyerr3=(randypred3-yt);
        % random forest tree
        rtree = fitrtree(x,y,'MaxNumSplits',100) %tunelength5, 100tree,
        ypred4 = predict(rtree,xt);
        yerr4=(ypred4-yt);
        randrtree = fitrtree(x,randy);
        randypred4 = predict(randrtree,xt);
        randyerr4=(randypred4-yt);

        %kernel
        % kernelmdl = fitrkernel(x,y,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
        % randMdl = fitrkernel(x,randy,'Learner','svm','NumExpansionDimensions',2048,'KernelScale',23.8858,'Lambda',2.9317e-04,'Epsilon',2.3795);
        % ypred4 = predict(kernelmdl,xt);
        % randypred4 = predict(randMdl,xt);
        % yerr4=mean(abs(ypred4-yt));
        % randyerr4=mean(abs(randypred4-yt));
        errmat=[errmat;[testid(:) yerr randyerr yerr2 randyerr2 yerr3 randyerr3 yerr4 randyerr4]]; %yerr4 randyerr4
        y_glm=[y_glm;ypred];
        y_las=[y_las;ypred2];
        y_svm=[y_svm;ypred3];
        y_rtree=[y_rtree;ypred4];
        % predictker=[predictker;[ypred4,yt]];
    end
    maes = mean(abs(errmat(:,2:end)));
    randcorr =[corr(errmat(:,2),errmat(:,3),'type','pearson'),corr(errmat(:,4),errmat(:,5),'type','pearson'),corr(errmat(:,6),errmat(:,7),'type','pearson'),corr(errmat(:,8),errmat(:,9),'type','pearson')];%,corr(errmat(:,8),errmat(:,9),'type','pearson')
    toc;
end