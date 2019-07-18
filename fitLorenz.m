%fitLorenz.m
%Edward Tekwa Nov 16, 2018
set(0,'defaultAxesFontSize',16)
TimeData=string(datetime);

%first, generate test data with observational noise
sigma = 10;
beta = 8/3;
rho = 28;
tEnd = 5000;
initXYZ=[1 1 1];
tRec=ceil(tEnd*0.25)+1; % time after which to record random points
numPts = 40; %number of random points to record per time-series
numExtrema = 20;
%numAvgPts = 100; %number of random points to use for averaging at the end of time-series
XYZindex=1; %set to vary and plot initial X, Y, or Z (1, 2, or 3)
%quant=numPts/(tEnd-tRec+1);
labels=['X','Y','Z'];

rhoRange=[0:.5:45]; %coarse grained
XinitRange=[-60:1:60]; %coarse grained
XYZOriginals=zeros(length(XinitRange),length(rhoRange),numPts,3);
XYZs=zeros(length(XinitRange),length(rhoRange),numPts,3);
XYZAvgs=zeros(length(XinitRange),length(rhoRange),3);
XYZExtrema=zeros(length(XinitRange),length(rhoRange),numExtrema,3);

noise=1; %measurement noise is additive in X, Y, and Z
XinitNoise=1; %noise in initial X, Y, and Z measurements
RhoNoise=1; %measurement noise in rho
XinitNoiseSteps=XinitNoise/(XinitRange(2)-XinitRange(1));
RhoNoiseSteps=RhoNoise/(rhoRange(2)-rhoRange(1));

for i=1:length(rhoRange)
    rho=rhoRange(i)
    for j=1:length(XinitRange)
        initXYZ(1)=XinitRange(j);
        f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
        [t,a] = ode45(f,[0 tEnd],initXYZ);     % Runge-Kutta 4th/5th order ODE solver
        [d,ix]=min(abs(t-[tRec:tEnd])); %find indices (ix) in time series corresponding to integer time steps
        States=a(datasample(ix,numPts,'Replace',false),:);
        XYZOriginals(j,i,:,:)=States; %record states near integer times
        XYZs(j,i,:,:)=States+noise*randn([numPts,3]); %add noise
        XYZtemp=a(ix,:);
        XYZAvgs(j,i,:)=mean(XYZtemp);
        diffXYZtemp=diff(XYZtemp);
        sortedDiffXYZtemp=sort(sum(abs(diffXYZtemp),2));
        candidateXYZs=XYZtemp(find(sum(abs(diffXYZtemp),2)<=sortedDiffXYZtemp(numExtrema+1))+1,:);
        XYZExtrema(j,i,:,:)=candidateXYZs(1:numExtrema,:);
    end
end

%Second, select random datapoints for empirical observational dataset
XYZindex=1;
sampleN=50; %number of independent time-series data
numMeanObs=2; %number of random unique observations to average out within time-series
numManyMeanObs=4;
%sampleRhorange=[1 length(rhoRange)]; %min and max observed rho (as indices in rhoRange)
%sampleRhorange=[1 find(rhoRange==25)-1]; %min and max observed r (as indices in rhoRange)
sampleRhorange=[find(rhoRange==25) length(rhoRange)]; %min and max observed r (as indices in rhoRange)
sampleRhoIndex=randi(sampleRhorange,1,sampleN); %array of size sampleN containing index of parameter r
sampleXinitIndex=randi([1 length(XinitRange)],1,sampleN); %array of size sampleN containing index of initial Y conditions
sampleRhoIndexWNoise=min(max(round(sampleRhoIndex+randn(size(sampleRhoIndex))*RhoNoiseSteps),1),length(rhoRange));
sampleXinitIndexWNoise=min(max(round(sampleXinitIndex+randn(size(sampleRhoIndex))*XinitNoiseSteps),1),length(XinitRange));
%sampleLag=randi([1 ceil(totTimeSteps*Period)],1,sampleN); %array of size sampleN containing time lag betwee initial and final observation
errS=zeros(1,sampleN); %record error^2 of each observation compared to prediction
errSmult=zeros(1,sampleN); %record error^2 of each observation compared to prediction
errSmanymult=zeros(1,sampleN); %record error^2 of each observation compared to prediction
Xobs=zeros(1,sampleN); %record Y observations
XmeanObs=zeros(1,sampleN); %record Ymean observations
XmanymeanObs=zeros(1,sampleN); %record Ymanymean observations
for i=1:sampleN
    samples=randperm(numPts,numManyMeanObs)'; %pick 'numManyMeanObs' random observations
    Xobs(i)=XYZs(sampleXinitIndexWNoise(i),sampleRhoIndexWNoise(i),samples(1),XYZindex); %only add noise here to distort observations
    XmeanObs(i)=mean(XYZs(sampleXinitIndexWNoise(i),sampleRhoIndexWNoise(i),samples(1:2),XYZindex)); %only add noise here to distort observations
    XmanymeanObs(i)=mean(XYZs(sampleXinitIndexWNoise(i),sampleRhoIndexWNoise(i),samples,XYZindex)); %only add noise here to distort observations
    errS(i)=(XYZAvgs(sampleXinitIndex(i),sampleRhoIndex(i),XYZindex)-Xobs(i))^2; %single observation
    errSmult(i)=(XYZAvgs(sampleXinitIndex(i),sampleRhoIndex(i),XYZindex)-XmeanObs(i))^2; %mean of multiple observations
    errSmanymult(i)=(XYZAvgs(sampleXinitIndex(i),sampleRhoIndex(i),XYZindex)-XmanymeanObs(i))^2; %mean of multiple observations
end
SS=sum(errS); %model sum of squares
SSmult=sum(errSmult); %model sum of squares
SSmanymult=sum(errSmanymult); %model sum of squares
SSMean=sum((Xobs-mean(Xobs)).^2); %mean sum of squares
SSMultMean=sum((XmeanObs-mean(XmeanObs)).^2); %mean sum of squares
SSManyMultMean=sum((XmanymeanObs-mean(XmanymeanObs)).^2); %mean sum of squares
R2=1-SS/SSMean
R2mult=1-SSmult/SSMultMean
R2manymult=1-SSmanymult/SSManyMultMean
Regress1=fitlm(rhoRange(sampleRhoIndex),Xobs);
Regress2=fitlm(rhoRange(sampleRhoIndex),XmeanObs);
Regress3=fitlm(rhoRange(sampleRhoIndex),XmanymeanObs);
RegressR2=Regress1.Rsquared.Ordinary
RegressMeanR2=Regress2.Rsquared.Ordinary
RegressManyMeanR2=Regress3.Rsquared.Ordinary
R2diff1=R2-RegressR2
R2diff2=R2mult-RegressMeanR2
R2diff3=R2manymult-RegressManyMeanR2


%permutation tests for significance
numRand=100000;
for j=1:numRand
    randRIndex=sampleRhoIndex(randperm(sampleN)); %randomize rho values
    randYinitIndex=sampleXinitIndex(randperm(sampleN)); %independently randomize initial conditions
    for i=1:sampleN
        %first, fit model to randomized r data
        RandRErrS(i)=(XYZAvgs(sampleXinitIndex(i),randRIndex(i),XYZindex)-Xobs(i))^2; %single observation
        RandRErrSmult(i)=(XYZAvgs(sampleXinitIndex(i),randRIndex(i),XYZindex)-XmeanObs(i))^2; %mean of multiple observations
        RandRErrSmanymult(i)=(XYZAvgs(sampleXinitIndex(i),randRIndex(i),XYZindex)-XmanymeanObs(i))^2; %mean of multiple observations
        %then, fit model to randomized Yinit data
        RandYiErrS(i)=(XYZAvgs(randYinitIndex(i),sampleRhoIndex(i),XYZindex)-Xobs(i))^2; %single observation
        RandYiErrSmult(i)=(XYZAvgs(randYinitIndex(i),sampleRhoIndex(i),XYZindex)-XmeanObs(i))^2; %mean of multiple observations
        RandYiErrSmanymult(i)=(XYZAvgs(randYinitIndex(i),sampleRhoIndex(i),XYZindex)-XmanymeanObs(i))^2; %mean of multiple observations  
    end
    RandR2(j)=1-sum(RandRErrS)/SSMean;
    RandR2mult(j)=1-sum(RandRErrSmult)/SSMultMean;
    RandR2manymult(j)=1-sum(RandRErrSmanymult)/SSManyMultMean;
    RandYiR2(j)=1-sum(RandYiErrS)/SSMean;
    RandYiR2mult(j)=1-sum(RandYiErrSmult)/SSMultMean;
    RandYiR2manymult(j)=1-sum(RandYiErrSmanymult)/SSManyMultMean;
    %RandRegressR2(j)
end
p_r1=(sum(RandR2>=R2)+1)/(length(RandR2)+1) %significance of bifurcation parameter modified by initial conditions
p_r2=(sum(RandR2mult>=R2mult)+1)/(length(RandR2mult)+1)
p_r3=(sum(RandR2manymult>=R2manymult)+1)/(length(RandR2manymult)+1)
p_Yi1=(sum(RandYiR2>=R2)+1)/(length(RandYiR2)+1) %significance of initial conditions alone
p_Yi2=(sum(RandYiR2mult>=R2mult)+1)/(length(RandYiR2mult)+1)
p_Yi3=(sum(RandYiR2manymult>=R2manymult)+1)/(length(RandYiR2manymult)+1)

scrsz = get(0,'ScreenSize');
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)]);
subplot(2,3,1)
scatter(reshape(repmat(rhoRange,length(XinitRange),1),1,[]),reshape(XYZAvgs(:,:,XYZindex),1,[]),0.1,'.k');
hold on
scatter(rhoRange(sampleRhoIndex),Xobs,30,'or','filled');
title (['single observations, R^2=' num2str(R2,2)])
xlabel '\rho'
ylabel 'X'
ylim([-20 20])
xlim(rhoRange(sampleRhorange))

subplot(2,3,2)
scatter(reshape(repmat(rhoRange,length(XinitRange),1),1,[]),reshape(XYZAvgs(:,:,XYZindex),1,[]),0.1,'.k');
hold on
scatter(rhoRange(sampleRhoIndex),XmeanObs,30,'or','filled');
title (['means of ' num2str(numMeanObs) ' observations, R^2=' num2str(R2mult,2)])
xlabel '\rho'
ylabel 'X'
ylim([-20 20])
xlim(rhoRange(sampleRhorange))

subplot(2,3,3)
scatter(reshape(repmat(rhoRange,length(XinitRange),1),1,[]),reshape(XYZAvgs(:,:,XYZindex),1,[]),0.1,'.k');
hold on
scatter(rhoRange(sampleRhoIndex),XmanymeanObs,30,'or','filled');
title (['means of ' num2str(numManyMeanObs) ' observations, R^2=' num2str(R2manymult,2)])
xlabel '\rho'
ylabel 'X'
ylim([-20 20])
xlim(rhoRange(sampleRhorange))

subplot(2,3,4)
%his=histogram(RandR2,20,'FaceColor',[0.9 0.9 0.9]);
hold on
[f,xi]=ksdensity(RandR2);
plot(xi,f,'LineWidth',2,'Color','r');
line([R2 R2], [0 max(f)],'LineWidth',2,'Color','k')
title (['p<' num2str(p_r1,2)])
xlabel 'R^2'
ylabel 'probability density'

subplot(2,3,5)
%his=histogram(RandR2mult,20,'FaceColor',[0.9 0.9 0.9]);
hold on
[f,xi]=ksdensity(RandR2mult);
plot(xi,f,'LineWidth',2,'Color','r');
line([R2mult R2mult], [0 max(f)],'LineWidth',2,'Color','k')
title (['p<' num2str(p_r2,2)])
xlabel 'R^2'
ylabel 'probability density'

subplot(2,3,6)
%his=histogram(RandR2manymult,20,'FaceColor',[0.9 0.9 0.9]);
hold on
[f,xi]=ksdensity(RandR2manymult);
plot(xi,f,'LineWidth',2,'Color','r');
line([R2manymult R2manymult], [0 max(f)],'LineWidth',2,'Color','k')
title (['p<' num2str(p_r3,2)]) %'10^%d'
xlabel 'R^2'
ylabel 'probability density'

FoodWebFile=sprintf('fitLorenz_%s.mat', TimeData);
save(FoodWebFile);