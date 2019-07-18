%Mumby.v
%Edward Tekwa June 21, 2019
%Fit Mumby coral solutions to linear data to test for false negative rate

%first, load coral-macroalgal bistable data and linear fit
%load('Mumby_single.mat') %contains LinModel as the linear model fitted to a single bistable data set

Mumby_single_linearFalsePos %run script to get initial parameter estimates for bistable model to a noisy linear data set
freeparstart=BestParamEsts;

%options = optimset('MaxFunEvals',2000,'MaxIter',2000,'Display','on','TolFun',1e-4,'TolX',1e-4);
options = optimset('MaxFunEvals',200,'MaxIter',200,'Display','off','TolFun',1e-4,'TolX',1e-4);
plotMarkSize=36;
%TimeData=string(datetime);

% Int=LinModel.Coefficients.Estimate(1);
% Slope=LinModel.Coefficients.Estimate(2);
% IntSD=LinModel.Coefficients.SE(1)*LinModel.NumObservations^0.5;
% SlopeSD=LinModel.Coefficients.SE(2)*LinModel.NumObservations^0.5;
% FixedIntSD=0.1;
% FixedSlopeSD=0.1;

%vectors to store R2 given different noisy linear datasets
numRep=1000;
FalsePosBistableRho2=nan(numRep,1);
FalsePosBistableR2=nan(numRep,1);
FalsePosLinearR2=nan(numRep,1);
FalsePosBistableParamEsts=nan(numRep,2);
FalsePosLinearParamEsts=nan(numRep,2);

for rep=1:numRep
    noisyregCMs=zeros(size(CMs));
    sampleCinit=[];
    sampleMinit=[];
    
    %regenerate data with underlying linear single attractor
    for pt=1:length(sampleGIndexWNoise) %for each g value in sample
        %regCMs(pt,1)=min(max(Int+IntSD*randn+(Slope+SlopeSD*randn)*g(sampleGIndexWNoise(pt)),0),1); %coral cover
        trueregCMs(pt,1)=min(max(Int+Slope*g(sampleGIndexWNoise(pt)),0),1); %coral cover
        noisyregCMs(pt,1)=min(max(Int+FixedIntSD*randn+(Slope+FixedSlopeSD*randn)*g(sampleGIndexWNoise(pt)),0),1); %coral cover
        noisyregCMs(pt,2)=(1-noisyregCMs(pt,1))*unifrnd(0,1,1); %macroalgal cover starts at any portion left over from coral cover
        sampleCinit(pt)=min(max((noisyregCMs(pt,1)-trueregCMs(pt,1))*2+trueregCMs(pt,1),0),1);
        sampleMinit(pt)=(1-sampleCinit(pt))*unifrnd(0,1,1); %initial macroalgal cover
    end
    
    CMs=noisyregCMs;
    
    freeparmin=[1 0 0 0 2]; % r,d,a,v,phi
    freeparmax=[1 1 0 2 Inf];
    BestParamEsts=[0 0 0 0 0];
    
    fvalMin=Inf; %reset minimum sum negative log likelihood
    estPhi=(mean(CMs(:,1))/(1-mean(CMs(:,1))))/var(CMs(:,1))-1; %approximate precision estimated from the normal mean of all data (disregarding relationship with bifurcation parameter)
    meanparstart=[mean(CMs(:,1)) estPhi]; %inital guesses for mode and precision
    meanparmin=[0 2];
    meanparmax=[Inf Inf];
    [MeanParamEsts,fval_mode,exitflag,output] = fminsearchbnd(@(params) Mean_LL_Beta(CMs(:,1)',params),meanparstart,meanparmin,meanparmax,options); %fit a constant line with beta distribution to data
    NegLL_mode=fval_mode/sampleN; %this is the null sum of negative likelihood
    NullT_Precision(T)=MeanParamEsts(2);
    NullT_fval(T)=fval_mode;
    
    %disp(['Sum of negative log likelihood and rho^2 in current search at t=' num2str(tRecord(T)) ':'])
    %for j=1:numFits
        %freeparstart=freeparmax.*unifrnd(0,1,[1 4]); %random initial guess
        %freeparstart=[freeparmax(1:4).*unifrnd(0,1,[1 4]) estPhi]; %random initial guess
        [ParamEsts,fval,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(sampleGIndexWNoise),CMs(:,1)',sampleCinit,sampleMinit,tEnd,params),freeparstart,freeparmin,freeparmax,options);
        Rho2FittedModel=1-exp(2*fval/sampleN)/exp(2*NegLL_mode);
        %disp([fval Rho2FittedModel])
        if fval<fvalMin
            BestParamEsts=ParamEsts;
            fvalMin=fval;
        end
    %end
    BestParamEsts(3)=BestParamEsts(4); %set estimated a=v
    %BestParamEsts
    BestRho2=1-exp(2*fvalMin/sampleN)/exp(2*NegLL_mode);
    SSMean=sum((CMs(:,1)-mean(CMs(:,1))).^2); %mean's sum of squares
    SSModel=Mumby_SS(g(sampleGIndexWNoise),CMs(:,1)',sampleCinit,sampleMinit,tEnd,BestParamEsts(1:4)); %model's sum of squares
    BestR2=1-SSModel/SSMean;
    FalsePosLinModel=fitlm(g(sampleGIndexWNoise),CMs(:,1)'); %,'intercept',false)
    R2FalsePosLinModel=FalsePosLinModel.Rsquared.Ordinary;
    
    FalsePosBistableR2(rep)=BestR2;
    FalsePosBistableRho2(rep)=BestRho2;
    FalsePosLinearR2(rep)=R2FalsePosLinModel;
    FalsePosBistableParamEsts(rep,:)=BestParamEsts([2,3]);
    FalsePosLinearParamEsts(rep,:)=FalsePosLinModel.Coefficients.Estimate';
    
    %summary statistics:
    disp(['LinearR2=' num2str(R2FalsePosLinModel,2) ', ModelR2=' num2str(BestR2,2) ', ModelRho2=' num2str(BestRho2,2)])
    
    
end
numRandDone=sum(~isnan(FalsePosBistableR2));
p_FalsePosLinearData=(sum(FalsePosLinearR2(1:numRandDone)<=FalsePosBistableR2(1:numRandDone)))/(numRandDone) %false positive rate (1-sided)
%p_FalsePosRhoLinearData=(sum(FalsePosLinearR2(1:numRandDone)<=FalsePosBistableRho2(1:numRandDone))+1)/(numRandDone+1) %false positive rate (1-sided)

mean(FalsePosBistableParamEsts)
std(FalsePosBistableParamEsts)
mean(FalsePosBistableR2)
std(FalsePosBistableR2)

mean(FalsePosLinearParamEsts)
std(FalsePosLinearParamEsts)
mean(FalsePosLinearR2)
std(FalsePosLinearR2)