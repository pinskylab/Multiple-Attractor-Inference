%Mumby_multipleRuns_parallel.m
%Edward Tekwa Nov 20, 18
%Compute Mumby coral solutions, fit and evaluate model on replicate sets

delete(gcp('nocreate'))
mycluster=parcluster;
parpool(min(mycluster.NumWorkers,16)) %run on either the max number of clusters or the limit specified here
options = optimset('MaxFunEvals',400,'MaxIter',400,'Display','off','TolFun',1e-3,'TolX',1e-3);
warning ('off','all');

TimeData=string(datetime);

numFits=5; %number of random initial guesses to try for each model fit
numRep=6; %number of replications per time lag
sampleN=40; %number of data samples per simulated dataset
tEnd=64; %25
tRecord=[0 1 2 4 8 16 32 64]; %record model fits at different elapsed times since initial disturbance
numRand=100; %number of permutations for significance tests

params=[1 0.5 1.1 1.1];
r=params(1); %coral overgrowth on algal turf
d=params(2); %coral mortality
a=params(3); %macroalgae overgrowth rate on coral
v=params(4); %gamma, macroalgae overgrowth rate on algal turf

g=[0:.005:1.2];

%known solutions for coral cover
Sol1=0*g;
Sol2_C_M=[[-(d - r + ((a + r)*(a*r^2 - 2*a^2*d - r*(a^4 + 2*a^3*r - 4*g*a^3 + 2*a^2*d*v + a^2*r^2 - 8*g*a^2*r + 4*g*a^2*v + 2*a*d*r*v - 4*g*a*r^2 + 4*g*a*r*v + d^2*v^2).^(1/2) + a^2*r - 2*a*d*r + 2*a*d*v + d*r*v))/(2*(a^3 + 2*a^2*r - v*a^2 + a*r^2 - v*a*r)))/r]; [(a*r^2 - 2*a^2*d - r*(a^4 + 2*a^3*r - 4*g*a^3 + 2*a^2*d*v + a^2*r^2 - 8*g*a^2*r + 4*g*a^2*v + 2*a*d*r*v - 4*g*a*r^2 + 4*g*a*r*v + d^2*v^2).^(1/2) + a^2*r - 2*a*d*r + 2*a*d*v + d*r*v)/(2*(a^3 + 2*a^2*r - v*a^2 + a*r^2 - v*a*r))]];
Sol2=Sol2_C_M(1,:);
Sol2(imag(Sol2)~=0)=-1;
Sol2_M=Sol2_C_M(2,:);
Sol2_M(imag(Sol2_M)~=0)=NaN;
Sol3=0*g+1-d/r;
Sols=[Sol1;Sol2;Sol3];

CinitRange=[0.01:0.01:0.8];

for i=1:length(g)
    if Sol2(i)>Sol3(i)
        Sols(3,i)=NaN;
        Sols(2,i)=NaN;
    end
    if Sol2(i)<Sol1(i)
        Sols(1,i)=NaN;
        Sols(2,i)=NaN;
    end
end



%perform fits on runs with different initial conditions
RepModelT_ParamEsts=zeros(length(tRecord),length(params)+1,numRep); %record parameter estimates at different elapsed times
RepModelT_fvalMin=zeros(length(tRecord),1,numRep);
RepModelT_R2Rho2=zeros(length(tRecord),2,numRep);
RepTrueT_fvalMin=zeros(length(tRecord),1,numRep);
RepTrueT_R2Rho2=zeros(length(tRecord),2,numRep);
RepLinT_R2=zeros(length(tRecord),1,numRep);
RepNullT_Precision=zeros(length(tRecord),1,numRep);
RepNullT_fval=zeros(length(tRecord),1,numRep);
RepModel_pValues=zeros(length(tRecord),7,numRep);
RepRand_Rho2=zeros(length(tRecord),3,numRep);

parfor Rep=1:numRep
    
    %generate noisy data
    tPlot=100; %plot at these time points
    Grange=[0.45 1];
    Cnoise=0.0; %beta random
    Gnoise=0.0; %normal random
    Mnoise=0.0; %assume no noise in macroalgae cover reading
    CnoiseSteps=Cnoise/(CinitRange(2)-CinitRange(1));
    GnoiseSteps=Gnoise/(g(2)-g(1));
    [e1 gMin]=min(abs(g-Grange(1))); %0.4
    [e2 gMax]=min(abs(g-Grange(2))); %0.9
    sampleGrange=[gMin gMax]; %min and max observed g (as indices). Eg. end:length(g)
    sampleGIndex=randi(sampleGrange,1,sampleN); %array of size sampleN containing index of parameter g
    sampleCinitIndex=randi([1 length(CinitRange)],1,sampleN); %array of size sampleN containing index of initial Y conditions
    sampleMinit=(1-CinitRange(sampleCinitIndex)).*unifrnd(0,1,1,sampleN); %macroalgal cover starts at any portion left over from coral cover
    sampleGIndexWNoise=min(max(round(sampleGIndex+randn(size(sampleGIndex))*GnoiseSteps),1),length(g)); %"observed" g
    sampleCinitIndexWNoise=min(max(round(sampleCinitIndex+randn(size(sampleGIndex))*CnoiseSteps),1),length(CinitRange)); %"observed Cinit"
    CMs=zeros(sampleN,2); %final states of coral and macroalgae with noise
    CMsTrue=zeros(sampleN,2); %true final states of coral and macroalgae
    CinitCond=zeros(sampleN,1); %initial condition of coral with noise
    CinitCondTrue=zeros(sampleN,1); %initial condition of coral
    [c1 index1]=min(abs(Sol2-Sol3(1)));
    [c2 index2]=min(abs(Sol2-Sol1(1)));
    
    t=[];
    X=[];
    XRecord=zeros(length(tRecord),sampleN,2);
    for i=1:sampleN
        gi=g(sampleGIndex(i));
        fMumby = @(t,X) [X(1)*(r*(1-X(1)-X(2))-d-a*X(2)); X(2)*(a*X(1)-gi/(X(2)+(1-X(1)-X(2)))+v*(1-X(1)-X(2)))]; %X(1) is coral cover, X(2) is macroalgal cover
        [t,X] = ode23(fMumby,[0 tEnd],[CinitRange(sampleCinitIndex(i)) sampleMinit(i)]);     % Runge-Kutta 2nd/3rd order ODE solver
        [t_diff t_index] = min(abs(tRecord-t)); %find nearest times in the time-series to record times
        XRecord(:,i,:)=X(t_index,:); %record coral and macroalgal time series
        if g(sampleGIndexWNoise(i))<g(index1) %region with sol1 (0) only
            CinitCond(i)=1; %basin of attraction of sol1
        elseif g(sampleGIndexWNoise(i))>g(index2) %region with sol3 only
            CinitCond(i)=2; %basin of attraction of sol3 (the second stable solution)
        elseif abs(X(end,1)-Sol1(sampleGIndexWNoise(i)))<abs(X(end,1)-Sol3(sampleGIndexWNoise(i))) %bistable region, projected long-run state close to lower stable state, below unstable solution
            CinitCond(i)=1;
        else %bistable region, above unstable solution
            CinitCond(i)=2;
        end
        t=[];
        X=[];
    end

    %fit model
    freeparmin=[1 0 0 0 2]; % r,d,a,v,phi
    freeparmax=[1 1 0 2 Inf];
    BestParamEsts=[0 0 0 0 0];
    
    ModelT_ParamEsts=zeros(length(tRecord),length(BestParamEsts)); %record parameter estimates at different elapsed times
    ModelT_fvalMin=zeros(length(tRecord),1);
    ModelT_R2Rho2=zeros(length(tRecord),2);
    TrueT_fvalMin=zeros(length(tRecord),1);
    TrueT_R2Rho2=zeros(length(tRecord),2);
    LinT_R2=zeros(length(tRecord),1);
    NullT_Precision=zeros(length(tRecord),1);
    NullT_fval=zeros(length(tRecord),1);
    Model_pValues=zeros(length(tRecord),8);
    Rand_Rho2=zeros(length(tRecord),3);
    
    for T=1:length(tRecord) %fit model at each time point in tRecord and store results
        %disp(['Rep=' num2str(Rep) ', t=' num2str(tRecord(T)) ':']);
        fvalMin=Inf; %reset minimum sum negative log likelihood
        CMs=reshape(XRecord(T,:,:),[],2); %record data nearest current elapsed time
        %fit and get negative log likelihood of the mean model (under beta distribution)
        estPhi=(mean(CMs(:,1))/(1-mean(CMs(:,1))))/var(CMs(:,1))-1; %approximate precision estimated from the normal mean of all data (disregarding relationship with bifurcation parameter)
        meanparstart=[mean(CMs(:,1)) estPhi]; %inital guesses for mode and precision
        meanparmin=[0 2];
        meanparmax=[Inf Inf];
        [MeanParamEsts,fval_mode,exitflag,output] = fminsearchbnd(@(params) Mean_LL_Beta(CMs(:,1)',params),meanparstart,meanparmin,meanparmax,options); %fit a constant line with beta distribution to data
        NegLL_mode=fval_mode/sampleN; %this is the null sum of negative likelihood
        NullT_Precision(T)=MeanParamEsts(2);
        NullT_fval(T)=fval_mode;
        
        %disp(['Sum of negative log likelihood and rho^2 in current search at t=' num2str(tRecord(T)) ':'])
        for j=1:numFits
            %freeparstart=freeparmax.*unifrnd(0,1,[1 4]); %random initial guess
            freeparstart=[freeparmax(1:4).*unifrnd(0,1,[1 4]) estPhi]; %random initial guess
            [ParamEsts,fval,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(sampleGIndexWNoise),CMs(:,1)',CinitRange(sampleCinitIndexWNoise),sampleMinit,tEnd,params),freeparstart,freeparmin,freeparmax,options);
            Rho2FittedModel=1-exp(2*fval/sampleN)/exp(2*NegLL_mode);
            %disp([fval Rho2FittedModel]);
            if fval<fvalMin
                BestParamEsts=ParamEsts;
                fvalMin=fval;
            end
        end
        BestParamEsts(3)=BestParamEsts(4); %set estimated a=v
        BestParamEsts;
        BestRho2=1-exp(2*fvalMin/sampleN)/exp(2*NegLL_mode);
        SSMean=sum((CMs(:,1)-mean(CMs(:,1))).^2); %mean's sum of squares
        SSModel=Mumby_SS(g(sampleGIndexWNoise),CMs(:,1)',CinitRange(sampleCinitIndexWNoise),sampleMinit,tEnd,BestParamEsts(1:4)); %model's sum of squares
        BestR2=1-SSModel/SSMean;
        LinModel=fitlm(g(sampleGIndexWNoise),CMs(:,1)'); %,'intercept',false)
        R2LinModel=LinModel.Rsquared.Ordinary;
        pLinModel=LinModel.Coefficients.pValue(2);
    
        %fix parameter values to true but estimate precision
        trueparstart=[params estPhi];
        trueparmin=[params 2];
        trueparmax=[params Inf];
        [TrueParamEsts,TrueFval,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(sampleGIndexWNoise),CMs(:,1)',CinitRange(sampleCinitIndexWNoise),sampleMinit,tEnd,params),trueparstart,trueparmin,trueparmax,options);
        TrueRho2=1-exp(2*TrueFval/sampleN)/exp(2*NegLL_mode);
        SSTrue=Mumby_SS(g(sampleGIndexWNoise),CMs(:,1)',CinitRange(sampleCinitIndexWNoise),sampleMinit,tEnd,params); %model's sum of squares
        TrueR2=1-SSTrue/SSMean;
        disp(['Rep=' num2str(Rep) ', t=' num2str(tRecord(T)) ': LinearR2=' num2str(R2LinModel,2) ', TrueR2=' num2str(TrueR2,2) ', TrueRho2=' num2str(TrueRho2,2) ', ModelR2=' num2str(BestR2,2) ', ModelRho2=' num2str(BestRho2,2)])
        disp(BestParamEsts)
        
        %permutation tests for significance
        RandGNegLL=nan(1,numRand);
        RandInitNegLL=nan(1,numRand);
        RandParamEsts1=nan(length(BestParamEsts),numRand);
        RandParamEsts2=nan(length(BestParamEsts),numRand);
        RandGRho2=nan(1,numRand);
        RandInitRho2=nan(1,numRand);
        RandGLinR2=nan(1,numRand);
        
        for j=1:numRand
            randGIndexWNoise=sampleGIndexWNoise(randperm(sampleN)); %randomize g values
            randInitIndex=randperm(sampleN);
            randCinitIndexWNoise=sampleCinitIndexWNoise(randInitIndex); %independently randomize initial conditions
            %first, fit model to randomized g data
            fval1=Inf;
            fval2=Inf;
            for k=1:numFits %search for parameters with this many different initial guesses
                if k==1 %start with best parameters for original data
                    tempparstart=BestParamEsts;
                else
                    tempparstart=[freeparmax(1:4).*unifrnd(0,1,[1 4]) estPhi]; %random initial guess
                end
                [tempParamEsts1,tempfval1,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(randGIndexWNoise),CMs(:,1)',CinitRange(sampleCinitIndexWNoise),sampleMinit,tEnd,params),tempparstart,freeparmin,freeparmax,options);
                %then, fit model to randomized Cinit data
                [tempParamEsts2,tempfval2,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(sampleGIndexWNoise),CMs(:,1)',CinitRange(randCinitIndexWNoise),sampleMinit(randInitIndex),tEnd,params),tempparstart,freeparmin,freeparmax,options);
                if tempfval1<fval1
                    ParamEsts1=tempParamEsts1;
                    fval1=tempfval1;
                end
                if tempfval2<fval2
                    ParamEsts2=tempParamEsts2;
                    fval2=tempfval2;
                end
                tempModelRho2=1-exp(2*fval1/sampleN)/exp(2*NegLL_mode);
                tempModel2Rho2=1-exp(2*fval2/sampleN)/exp(2*NegLL_mode);
                tempLinModel=fitlm(g(randGIndexWNoise),CMs(:,1)'); %,'intercept',false)
                tempR2LinModel=tempLinModel.Rsquared.Ordinary;
            end
            RandGNegLL(j)=fval1;
            RandInitNegLL(j)=fval2;
            RandParamEsts1(:,j)=ParamEsts1';
            RandParamEsts2(:,j)=ParamEsts2';
            RandGRho2(j)=tempModelRho2;
            RandInitRho2(j)=tempModel2Rho2;
            RandGLinR2(j)=tempR2LinModel;
            disp(['Rep=' num2str(Rep) ', t=' num2str(tRecord(T)) ': rand' num2str(j) ', RandLinearR2=' num2str(tempR2LinModel,2) ', RandGRho2=' num2str(tempModelRho2,2) ', RandInitRho2=' num2str(tempModel2Rho2,2)])
        end
        numRandDone=sum(~isnan(RandGNegLL));
        p_G=(sum(RandGNegLL(1:numRandDone)<=fvalMin)+1)/(numRandDone+1); %significance of bifurcation parameter modified by initial conditions (1-sided)
        p_IC=(sum(RandInitNegLL(1:numRandDone)<=fvalMin)+1)/(numRandDone+1); %significance of initial conditions (1-sided)
        p_Model_Lin=(sum(RandGRho2(1:numRandDone)-RandGLinR2(1:numRandDone)>=BestRho2-R2LinModel)+1)/(numRandDone+1); %significance of model against linear regression (1-sided)
        p_d_G=(2*min(sum(RandParamEsts1(2,1:numRandDone)>=BestParamEsts(2)),sum(RandParamEsts1(2,1:numRandDone)<=BestParamEsts(2)))+1)/(numRandDone+1); %significance of parameter d aganist randomized bifurcation parameter g (2-sided)
        p_v_G=(2*min(sum(RandParamEsts1(4,1:numRandDone)>=BestParamEsts(4)),sum(RandParamEsts1(4,1:numRandDone)<=BestParamEsts(4)))+1)/(numRandDone+1); %significance of parameter v aganist randomized bifurcation parameter g (2-sided)
        p_d_IC=(2*min(sum(RandParamEsts2(2,1:numRandDone)>=BestParamEsts(2)),sum(RandParamEsts2(2,1:numRandDone)<=BestParamEsts(2)))+1)/(numRandDone+1); %significance of parameter d aganist randomized initial condition (2-sided)
        p_v_IC=(2*min(sum(RandParamEsts2(4,1:numRandDone)>=BestParamEsts(4)),sum(RandParamEsts2(4,1:numRandDone)<=BestParamEsts(4)))+1)/(numRandDone+1); %significance of parameter v aganist randomized initial condition (2-sided)
        mean_RandGLinR2=mean(RandGLinR2);
        mean_RandGRho2=mean(RandGRho2);
        mean_RandInitRho2=mean(RandInitRho2);
        
        %summary statistics:
        disp(['Rep=' num2str(Rep) ', t=' num2str(tRecord(T)) ': pG=' num2str(p_G,2) ' ,pIC=' num2str(p_IC,2) ' ,p_Model_Lin=' num2str(p_Model_Lin,2) ])
        %disp(['Parameter estimates: ' BestParamEsts])
        
        
        ModelT_ParamEsts(T,:)=BestParamEsts;
        ModelT_fvalMin(T)=fvalMin;
        ModelT_R2Rho2(T,:)=[BestR2 BestRho2];
        TrueT_fvalMin(T)=TrueFval;
        TrueT_R2Rho2(T,:)=[TrueR2 TrueRho2];
        LinT_R2(T)=R2LinModel;
        Model_pValues(T,:)=[p_G p_IC p_Model_Lin p_d_G p_v_G p_d_IC p_v_IC ];
        Rand_Rho2(T,:)=[mean_RandGRho2 mean_RandInitRho2 mean_RandGLinR2];
    end
    RepModelT_ParamEsts(:,:,Rep)=ModelT_ParamEsts;
    RepModelT_fvalMin(:,:,Rep)=ModelT_fvalMin;
    RepModelT_R2Rho2(:,:,Rep)=ModelT_R2Rho2;
    RepTrueT_fvalMin(:,:,Rep)=TrueT_fvalMin;
    RepTrueT_R2Rho2(:,:,Rep)=TrueT_R2Rho2;
    RepLinT_R2(:,:,Rep)=LinT_R2;
    RepNullT_Precision(:,:,Rep)=NullT_Precision;
    RepNullT_fval(:,:,Rep)=NullT_fval;
    RepModel_pValues(:,:,Rep)=Model_pValues;
    RepRand_Rho2(:,:,Rep)=Rand_Rho2;
end

%get mean and bootstrapped 95% confidence intervals for p-values of the
%model and of model compared to linear regression
mean_p_G=mean(RepModel_pValues(:,1,:),3)
mean_p_Model_Lin=mean(RepModel_pValues(:,3,:),3);
sd_p_G=std(RepModel_pValues(:,1,:),0,3)
sd_p_Model_Lin=std(RepModel_pValues(:,3,:),0,3);
up_p_G=zeros(size(RepModel_pValues,1),1);
lo_p_G=zeros(size(RepModel_pValues,1),1);
up_p_Model_Lin=zeros(size(RepModel_pValues,1),1);
lo_p_Model_Lin=zeros(size(RepModel_pValues,1),1);

falsPos010=cdf('normal',0.1,mean_p_G(1),sd_p_G(1))
falsPos005=cdf('normal',0.05,mean_p_G(1),sd_p_G(1))
falsPos001=cdf('normal',0.01,mean_p_G(1),sd_p_G(1))

MumbyFileName=sprintf('MumbyFit_%s.mat', TimeData);
save(MumbyFileName);

%plot model fit results at different elapsed times

set(0,'defaultAxesFontSize',16)
set(0,'DefaultLineLineWidth',1)
scrsz = get(0,'ScreenSize');
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.5 scrsz(4)/3]);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
subplot(1,2,1)
hold on
yyaxis left
bl1=boundedline([1:length(tRecord)],nanmean((reshape(RepModelT_ParamEsts(:,2,:),[],size(RepModelT_ParamEsts,3)))'),[nanstd((reshape(RepModelT_ParamEsts(:,2,:),[],size(RepModelT_ParamEsts,3)))')*1.96./((sum(~isnan((reshape(RepModelT_ParamEsts(:,2,:),[],size(RepModelT_ParamEsts,3)))'))).^0.5)],'b','alpha'); drawnow; set(bl1,'linewidth',2);
bl1=boundedline([1:length(tRecord)],nanmean((reshape(RepModelT_ParamEsts(:,4,:),[],size(RepModelT_ParamEsts,3)))'),[nanstd((reshape(RepModelT_ParamEsts(:,4,:),[],size(RepModelT_ParamEsts,3)))')*1.96./((sum(~isnan((reshape(RepModelT_ParamEsts(:,4,:),[],size(RepModelT_ParamEsts,3)))'))).^0.5)],'r','alpha'); drawnow; set(bl1,'linewidth',2);
%plot(reshape(RepModelT_ParamEsts(:,2,:),[],size(RepModelT_ParamEsts,3)),'-b'); %estimated d
%plot(reshape(RepModelT_ParamEsts(:,4,:),[],size(RepModelT_ParamEsts,3)),'-r'); %estimated v
ref1=refline(0,params(2));
ref1.LineStyle='--';
ref1.Color='b';
ref2=refline(0,params(4));
ref2.LineStyle='--';
ref2.Color='r';
ylabel 'parameters'
yyaxis right
bl1=boundedline([1:length(tRecord)],nanmean((reshape(log2(RepModelT_ParamEsts(:,5,:)),[],size(RepModelT_ParamEsts,3)))'),[nanstd((reshape(log2(RepModelT_ParamEsts(:,5,:)),[],size(RepModelT_ParamEsts,3)))')*1.96./((sum(~isnan((reshape(log2(RepModelT_ParamEsts(:,5,:)),[],size(RepModelT_ParamEsts,3)))'))).^0.5)],'k','alpha'); drawnow; set(bl1,'linewidth',2);
%plot(reshape(log2(RepModelT_ParamEsts(:,5,:)),[],size(RepModelT_ParamEsts,3)),'-k'); %estimated precision
ylabel 'log_2(precision)'
xlabel 'final observation time'
xticks([1:length(tRecord)])
xticklabels(tRecord)
xlim([1 length(tRecord)])

subplot(1,2,2)
hold on
yyaxis right
bl1=boundedline([1:length(tRecord)],nanmean((reshape(RepModelT_R2Rho2(:,1,:),[],size(RepModelT_ParamEsts,3)))'),[nanstd((reshape(RepModelT_R2Rho2(:,1,:),[],size(RepModelT_ParamEsts,3)))')*1.96./((sum(~isnan((reshape(RepModelT_R2Rho2(:,1,:),[],size(RepModelT_ParamEsts,3)))'))).^0.5)],'-r','alpha'); drawnow; set(bl1,'linewidth',2); %model R2
bl1=boundedline([1:length(tRecord)],nanmean((reshape(RepModelT_R2Rho2(:,2,:),[],size(RepModelT_ParamEsts,3)))'),[nanstd((reshape(RepModelT_R2Rho2(:,2,:),[],size(RepModelT_ParamEsts,3)))')*1.96./((sum(~isnan((reshape(RepModelT_R2Rho2(:,2,:),[],size(RepModelT_ParamEsts,3)))'))).^0.5)],'--r','alpha'); drawnow; set(bl1,'linewidth',2); %model Rho2
%plot(reshape(RepModelT_R2Rho2(:,1,:),[],size(RepModelT_ParamEsts,3)),'-r'); %model R2
%plot(reshape(RepModelT_R2Rho2(:,2,:),[],size(RepModelT_ParamEsts,3)),'--r'); %model Rho2
bl1=boundedline([1:length(tRecord)],nanmean((reshape(RepLinT_R2,[],size(RepModelT_ParamEsts,3)))'),[nanstd((reshape(RepLinT_R2,[],size(RepModelT_ParamEsts,3)))')*1.96./((sum(~isnan((reshape(RepLinT_R2,[],size(RepModelT_ParamEsts,3)))'))).^0.5)],'-b','alpha'); drawnow; set(bl1,'linewidth',2); %linear R2
%plot(reshape(RepLinT_R2,[],size(RepModelT_ParamEsts,3)),'-b'); %linear R2
xlabel 'final observation time'
ylabel 'fit'
xticks([1:length(tRecord)])
xticklabels(tRecord)
xlim([1 length(tRecord)])

yyaxis left
bl1=boundedline([1:length(tRecord)],nanmean((reshape(RepModel_pValues(:,1,:),[],size(RepModelT_ParamEsts,3)))'),[nanstd((reshape(RepModel_pValues(:,1,:),[],size(RepModelT_ParamEsts,3)))')*1.96./((sum(~isnan((reshape(RepModel_pValues(:,1,:),[],size(RepModelT_ParamEsts,3)))'))).^0.5)],'-k','alpha'); drawnow; set(bl1,'linewidth',2); %model p-value
%plot(reshape(RepModel_pValues(:,1,:),[],size(RepModel_pValues,3)),'-k'); %model p-value
ref1=refline(0,0.05);
ref1.LineStyle='-';
ref1.Color='m';
ylabel 'p-value'