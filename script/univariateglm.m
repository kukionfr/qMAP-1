function pmresuf = univariateglm(cfmat,y0,runtag,flbl,ccm)
    pmresuf=[];
    for kff=1:size(cfmat,2) %iterate each feature for univariate 
        x0=cfmat(:,kff);
        yerrs=[];ypreds=[];yts=[];
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
            yerr=abs(ypred-yt);
            yerrs =[yerrs;yerr];
            ypreds=[ypreds;ypred];
            yts =[yts,yt];
        end 
        pmresuf=[pmresuf;mean(yerrs) std(yerrs) mean(yerrs(ccm)) mean(yerrs(~ccm))];
    end
    xlsname_unipred=['../output/Univariate Predicting power',runtag,'.xlsx'];
    tmp0=[flbl(:,:); num2cell(pmresuf')]';
    tmp0lbl={'feature id','tissue area','feature group','variable name','error-mean','error-std','error-male','error-female'};
    tmp0=[tmp0lbl;tmp0];
    xlswrite(xlsname_unipred,tmp0);
end