%Mumby.v
%Edward Tekwa June 21, 2019
%Fit Mumby coral solutions to linear data to test for false negative rate

%first, load coral-macroalgal bistable data and linear fit
load('Mumby_single.mat') %contains LinModel as the linear model fitted to a single bistable data set

%options = optimset('MaxFunEvals',2000,'MaxIter',2000,'Display','on','TolFun',1e-4,'TolX',1e-4);
options = optimset('MaxFunEvals',200,'MaxIter',200,'Display','off','TolFun',1e-4,'TolX',1e-4);
plotMarkSize=36;
TimeData=string(datetime);

Int=LinModel.Coefficients.Estimate(1);
Slope=LinModel.Coefficients.Estimate(2);
IntSD=LinModel.Coefficients.SE(1)*LinModel.NumObservations^0.5;
SlopeSD=LinModel.Coefficients.SE(2)*LinModel.NumObservations^0.5;
FixedIntSD=0.1;
FixedSlopeSD=0.1;

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
    sampleMinit(pt)=(1-sampleCinit(pt))*unifrnd(0,1,1);
    %macroalgal cover
end

CMs=noisyregCMs;

figure;
scatter(g(sampleGIndexWNoise),CMs(:,1),'k');
hold on
scatter(g(sampleGIndexWNoise),sampleCinit,'r');

%fit model
%freeparmin=[1 0 1 0 2]; % r,d,a,v,phi
%freeparmax=[1 1 1 2 Inf];
freeparmin=[1 0 0 0 2]; % r,d,a,v,phi
freeparmax=[1 1 0 2 Inf];
%fvalMin=Inf;
BestParamEsts=[0 0 0 0 0];

% ModelT_ParamEsts=zeros(length(tRecord),length(BestParamEsts)); %record parameter estimates at different elapsed times
% ModelT_fvalMin=zeros(length(tRecord),1);
% ModelT_R2Rho2=zeros(length(tRecord),2);
% TrueT_fvalMin=zeros(length(tRecord),1);
% TrueT_R2Rho2=zeros(length(tRecord),2);
% LinT_R2=zeros(length(tRecord),1);
% NullT_Precision=zeros(length(tRecord),1);
% NullT_fval=zeros(length(tRecord),1);

%for T=1:length(tRecord) %fit model at each time point in tRecord and store results
fvalMin=Inf; %reset minimum sum negative log likelihood
%CMsTrue=reshape(XRecord(T,:,:),[],2); %record data nearest current elapsed time
%CMs=max(reshape(XRecord(T,:,:),[],2)+Cnoise*(randn(1,2)),0);
%CMs=reshape(XRecord(T,:,:),[],2); %record data nearest current elapsed time
%fit and get negative log likelihood of the mean model (under beta distribution)
estPhi=(mean(CMs(:,1))/(1-mean(CMs(:,1))))/var(CMs(:,1))-1; %approximate precision estimated from the normal mean of all data (disregarding relationship with bifurcation parameter)
meanparstart=[mean(CMs(:,1)) estPhi]; %inital guesses for mode and precision
meanparmin=[0 2];
meanparmax=[Inf Inf];
[MeanParamEsts,fval_mode,exitflag,output] = fminsearchbnd(@(params) Mean_LL_Beta(CMs(:,1)',params),meanparstart,meanparmin,meanparmax,options); %fit a constant line with beta distribution to data
NegLL_mode=fval_mode/sampleN; %this is the null sum of negative likelihood
NullT_Precision(T)=MeanParamEsts(2);
NullT_fval(T)=fval_mode;

disp(['Sum of negative log likelihood and rho^2 in current search at t=' num2str(tRecord(T)) ':'])
for j=1:numFits
    %freeparstart=freeparmax.*unifrnd(0,1,[1 4]); %random initial guess
    freeparstart=[freeparmax(1:4).*unifrnd(0,1,[1 4]) estPhi]; %random initial guess
    [ParamEsts,fval,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(sampleGIndexWNoise),CMs(:,1)',sampleCinit,sampleMinit,tEnd,params),freeparstart,freeparmin,freeparmax,options);
    Rho2FittedModel=1-exp(2*fval/sampleN)/exp(2*NegLL_mode);
    disp([fval Rho2FittedModel])
    if fval<fvalMin
        BestParamEsts=ParamEsts;
        fvalMin=fval;
    end
end
BestParamEsts(3)=BestParamEsts(4); %set estimated a=v
BestParamEsts
BestRho2=1-exp(2*fvalMin/sampleN)/exp(2*NegLL_mode);
SSMean=sum((CMs(:,1)-mean(CMs(:,1))).^2); %mean's sum of squares
SSModel=Mumby_SS(g(sampleGIndexWNoise),CMs(:,1)',sampleCinit,sampleMinit,tEnd,BestParamEsts(1:4)); %model's sum of squares
BestR2=1-SSModel/SSMean;
FalsePosLinModel=fitlm(g(sampleGIndexWNoise),CMs(:,1)') %,'intercept',false)
R2FalsePosLinModel=FalsePosLinModel.Rsquared.Ordinary;


%True_LL=Mumby_LL(g(sampleGIndexWNoise),CMs(:,1)',sampleCinit,params); %check LL for true model
%TrueR2=1-True_LL/SSMean
%fix parameter values to true but estimate precision
trueparstart=[params estPhi];
trueparmin=[params 2];
trueparmax=[params Inf];
[TrueParamEsts,TrueFval,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(sampleGIndexWNoise),CMs(:,1)',sampleCinit,sampleMinit,tEnd,params),trueparstart,trueparmin,trueparmax,options);
TrueRho2=1-exp(2*TrueFval/sampleN)/exp(2*NegLL_mode);
SSTrue=Mumby_SS(g(sampleGIndexWNoise),CMs(:,1)',sampleCinit,sampleMinit,tEnd,params); %model's sum of squares
TrueR2=1-SSTrue/SSMean;

%summary statistics:
disp(['LinearRho2=' num2str(R2FalsePosLinModel,2) ', TrueR2=' num2str(TrueR2,2) ', TrueRho2=' num2str(TrueRho2,2) ', ModelR2=' num2str(BestR2,2) ', ModelRho2=' num2str(BestRho2,2)])

%if sum(tRecord(T)==tPlot)==1 %plot model
set(0,'defaultAxesFontSize',16)
set(0,'defaultaxeslinewidth',2)
scrsz = get(0,'ScreenSize');
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1.5 scrsz(4)/2.5]);
subplot(1,2,1)
hold on
cmap=colormap(parula(128)); % set colormap
color_i = fix(sampleMinit*127)+1;
Minit_colors=cmap(color_i,:);
lowInitIndices=(CinitCond==1);
hiInitIndices=(CinitCond==2);
grazing=g(sampleGIndex);
Cinit=sampleCinit;
%gscatter(g(sampleGIndexWNoise),CMs(:,1),CinitCond,'kkk','x.o'); %plot final states with true initial conditions
%gscatter(g(sampleGIndex),CinitRange(sampleCinitIndex),CinitCond,'kkk','o.x'); %plot initial states with true initial conditions
%gscatter(g(sampleGIndex),CinitRange(sampleCinitIndex),CinitCond,Minit_colors,'o.x'); %plot initial states with true initial conditions
scatter(grazing(lowInitIndices),Cinit(lowInitIndices),plotMarkSize,Minit_colors(lowInitIndices,:)); %plot final states with modelled initial conditions
scatter(grazing(hiInitIndices),Cinit(hiInitIndices),plotMarkSize,Minit_colors(hiInitIndices,:),'filled'); %plot final states with modelled initial conditions

%sampleMinit
plot(g,Sols(1,:),'k','lineWidth',2);
hold on
plot(g,Sols(2,:),'--k','lineWidth',2);
plot(g,Sols(3,:),'k','lineWidth',2);
%ylim([-0.1 1])
ylim([min(CinitRange)-0.05 max(CinitRange)+0.05])
xlim([Grange(1)-0.05 Grange(2)+0.05])
xlabel 'grazing rate'
ylabel 'initial coral cover'
title ''

subplot(1,2,2)
hold on
qplot=quiver(g(sampleGIndex),sampleCinit,zeros(1,sampleN),CMs(:,1)'-sampleCinit,0); %plot arrows from initial to final states
qplot.Color='black';
qplot.MaxHeadSize=0.1;
%legend 'C_0\epsilon no coral' 'C_0\epsilon high coral'

r=BestParamEsts(1); %coral overgrowth on algal turf
d=BestParamEsts(2); %coral mortality
a=BestParamEsts(3); %macroalgae overgrowth rate on coral
v=BestParamEsts(4); %gamma, macroalgae overgrowth rate on algal turf
phi=BestParamEsts(5); %phi, precision of beta distributed coral cover

%known solutions for coral cover
Sol1=0*g;
%Sol2_1=-(d-r+(r+v)*(r*v-d*v+v^2+((v*(v*d^2+2*d*r*v+2*d*v^2+v*r^2-4*g*r^2+2*r*v^2-4*g*r*v+v^3)).^(.5)))/(2*(v^2+r*v)))/r;
Sol2_C_M=[[-(d - r + ((a + r)*(a*r^2 - 2*a^2*d - r*(a^4 + 2*a^3*r - 4*g*a^3 + 2*a^2*d*v + a^2*r^2 - 8*g*a^2*r + 4*g*a^2*v + 2*a*d*r*v - 4*g*a*r^2 + 4*g*a*r*v + d^2*v^2).^(1/2) + a^2*r - 2*a*d*r + 2*a*d*v + d*r*v))/(2*(a^3 + 2*a^2*r - v*a^2 + a*r^2 - v*a*r)))/r]; [(a*r^2 - 2*a^2*d - r*(a^4 + 2*a^3*r - 4*g*a^3 + 2*a^2*d*v + a^2*r^2 - 8*g*a^2*r + 4*g*a^2*v + 2*a*d*r*v - 4*g*a*r^2 + 4*g*a*r*v + d^2*v^2).^(1/2) + a^2*r - 2*a*d*r + 2*a*d*v + d*r*v)/(2*(a^3 + 2*a^2*r - v*a^2 + a*r^2 - v*a*r))]];
Sol2=Sol2_C_M(1,:);
Sol2(imag(Sol2)~=0)=-1;
Sol2_M=Sol2_C_M(2,:);
Sol2_M(imag(Sol2_M)~=0)=NaN;
Sol3=0*g+1-d/r;
Sols=[Sol1;Sol2;Sol3];

%CinitRange=[0:.01:Sol3(1)*1.2];

%         for i=1:length(g)
%             if Sol2(i)>Sol3(i)
%                 Sols(3,i)=NaN;
%                 Sols(2,i)=NaN;
%             end
%             if Sol2(i)<Sol1(i)
%                 Sols(1,i)=NaN;
%                 Sols(2,i)=NaN;
%             end
%         end
%         [c1 index1]=min(abs(Sol2-Sol3(1)));
%         [c2 index2]=min(abs(Sol2-Sol1(1)));

for i=1:sampleN
    %simulate long-term outcomes using estimated parameters to decide
    %basin of attraction
    fMumby = @(t,X) [X(1)*(r*(1-X(1)-X(2))-d-a*X(2)); X(2)*(a*X(1)-g(sampleGIndexWNoise(i))/(X(2)+(1-X(1)-X(2)))+v*(1-X(1)-X(2)))]; %X(1) is coral cover, X(2) is macroalgal cover
    [t,X] = ode45(fMumby,[0 tEnd],[CinitRange(sampleCinitIndex(i)) sampleMinit(i)]);     % Runge-Kutta 4th/5th order ODE solver
    
    if Sol2(sampleGIndexWNoise(i))>Sol3(sampleGIndexWNoise(i)) %solution is 0
        CinitCond(i)=1; %basin of attraction of sol1
    elseif Sol2(sampleGIndexWNoise(i))<Sol1(sampleGIndexWNoise(i)) %solution is the higher equilibrium
        CinitCond(i)=2; %basin of attraction of sol3 (the second stable solution)
        %elseif CinitRange(sampleCinitIndexWNoise(i))/sampleMinit(i)<Sol2(sampleGIndexWNoise(i))/Sol2_M(sampleGIndexWNoise(i))  %bistable region, below unstable solution
    elseif abs(X(end,1)-Sol1(sampleGIndexWNoise(i)))<abs(X(end,1)-Sol3(sampleGIndexWNoise(i)))  %bistable region, projected long-run state close to lower stable state, below unstable solution
        CinitCond(i)=1;
        %elseif CinitRange(sampleCinitIndexWNoise(i))/sampleMinit(i)==Sol2(sampleGIndexWNoise(i))/Sol2_M(sampleGIndexWNoise(i))  %on unstable solution
        %    CinitCond(i)=3;
    else %bistable region, above unstable solution
        CinitCond(i)=2;
    end
end
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

%fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/2]);
plot(g,Sols(1,:),'r','lineWidth',2);
hold on
plot(g,Sols(2,:),'--r','lineWidth',2);
plot(g,Sols(3,:),'r','lineWidth',2);
%gscatter(g(sampleGIndex),CMs(:,1),CinitCond,'rrr','o.x'); %plot final states with modelled initial conditions
color_i = fix(CMs(:,2)*127)+1;
Mfinal_colors=cmap(color_i,:);
lowInitIndices=(CinitCond==1);
hiInitIndices=(CinitCond==2);
grazing=g(sampleGIndex);
scatter(grazing(lowInitIndices),CMs(lowInitIndices,1),plotMarkSize,Mfinal_colors(lowInitIndices,:)); %plot final states with modelled initial conditions
scatter(grazing(hiInitIndices),CMs(hiInitIndices,1),plotMarkSize,Mfinal_colors(hiInitIndices,:),'filled'); %plot final states with modelled initial conditions
ylim([min(CinitRange)-0.05 max(CinitRange)+0.05])
xlim([Grange(1)-0.05 Grange(2)+0.05])

LinPlot=plot(FalsePosLinModel);
set(LinPlot,'Marker','none','LineWidth',2,'Color','k')
legend off
xlabel 'grazing rate'
ylabel 'final coral cover'
title(['N=' num2str(sampleN) ', tRun=' num2str(tRecord(T)) ', linR2=' num2str(R2FalsePosLinModel,2) ', tR2=' num2str(TrueR2,2) ', tRho2=' num2str(TrueRho2,2) ', mR2=' num2str(BestR2,2) ', mRho2=' num2str(BestRho2,2)])
%end

%     ModelT_ParamEsts(T,:)=BestParamEsts;
%     ModelT_fvalMin(T)=fvalMin;
%     ModelT_R2Rho2(T,:)=[BestR2 BestRho2];
%     TrueT_fvalMin(T)=TrueFval;
%     TrueT_R2Rho2(T,:)=[TrueR2 TrueRho2];
%     LinT_R2(T)=R2FalsePosLinModel;
%end

% %plot model fit results at different elapsed times
% fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/2]);
% subplot(1,2,1)
% hold on
% yyaxis left
% plot(ModelT_ParamEsts(:,2),'b'); %estimated d
% plot(ModelT_ParamEsts(:,4),'r'); %estimated v
% ylabel 'parameters'
% yyaxis right
% plot(log2(ModelT_ParamEsts(:,5)),'k'); %estimated precision
% plot(log2(NullT_Precision),'g'); %estimated precision
% ylabel 'log_2(precision)'
% xlabel 'final observation time'
%
% subplot(1,2,2)
% hold on
% plot(ModelT_R2Rho2(:,1),'k'); %model R2
% plot(ModelT_R2Rho2(:,2),'r'); %model Rho2
% plot(LinT_R2,'b'); %model Rho2
% xlabel 'final observation time'
% ylabel 'fit'

%
% %permutation tests for significance
% numRand=100;
% RandGNegLL=nan(1,numRand);
% RandInitNegLL=nan(1,numRand);
% RandParamEsts1=nan(length(BestParamEsts),numRand);
% RandParamEsts2=nan(length(BestParamEsts),numRand);
% RandGRho2=nan(1,numRand);
% RandInitRho2=nan(1,numRand);
% RandGLinR2=nan(1,numRand);
%
% for j=1:numRand
%     randGIndexWNoise=sampleGIndexWNoise(randperm(sampleN)); %randomize g values
%     randInitIndex=randperm(sampleN);
%     randCinitIndexWNoise=sampleCinitIndexWNoise(randInitIndex); %independently randomize initial conditions
%     %first, fit model to randomized g data
%     %[ParamEsts,fval1,exitflag,output] = fminsearchbnd(@(params) Mumby_LL(g(randGIndexWNoise),CMs(:,1)',sampleCinit,params),BestParamEsts,freeparmin,freeparmax,options);
%     %[ParamEsts1,fval1,exitflag,output] = fminsearchbnd(@(params) Mumby_LL(g(randGIndexWNoise),CMs(:,1)',sampleCinit,params),BestParamEsts,freeparmin,freeparmax,options);
%     fval1=Inf;
%     fval2=Inf;
%     for k=1:numFits %search for parameters with this many different initial guesses
%         if k==1 %start with best parameters for original data
%             tempparstart=BestParamEsts;
%         else
%             tempparstart=[freeparmax(1:4).*unifrnd(0,1,[1 4]) estPhi]; %random initial guess
%         end
%         [tempParamEsts1,tempfval1,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(randGIndexWNoise),CMs(:,1)',sampleCinit,sampleMinit,tEnd,params),tempparstart,freeparmin,freeparmax,options);
%         %then, fit model to randomized Cinit data
%         %[ParamEsts2,fval2,exitflag,output] = fminsearchbnd(@(params) Mumby_LL(g(sampleGIndexWNoise),CMs(:,1)',CinitRange(randCinitIndexWNoise),params),BestParamEsts,freeparmin,freeparmax,options);
%         [tempParamEsts2,tempfval2,exitflag,output] = fminsearchbnd(@(params) Mumby_LL_Beta(g(sampleGIndexWNoise),CMs(:,1)',CinitRange(randCinitIndexWNoise),sampleMinit(randInitIndex),tEnd,params),tempparstart,freeparmin,freeparmax,options);
%         if tempfval1<fval1
%             ParamEsts1=tempParamEsts1;
%             fval1=tempfval1;
%         end
%         if tempfval2<fval2
%             ParamEsts2=tempParamEsts2;
%             fval2=tempfval2;
%         end
%         tempModelRho2=1-exp(2*fval1/sampleN)/exp(2*NegLL_mode);
%         tempModel2Rho2=1-exp(2*fval2/sampleN)/exp(2*NegLL_mode);
%         tempFalsePosLinModel=fitlm(g(randGIndexWNoise),CMs(:,1)'); %,'intercept',false)
%         tempR2FalsePosLinModel=tempFalsePosLinModel.Rsquared.Ordinary;
%     end
%     disp(['rand' num2str(j) ', RandLinearR2=' num2str(tempR2FalsePosLinModel,2) ', RandGRho2=' num2str(tempModelRho2,2) ', RandInitRho2=' num2str(tempModel2Rho2,2)])
%     RandGNegLL(j)=fval1;
%     RandInitNegLL(j)=fval2;
%     RandParamEsts1(:,j)=ParamEsts1';
%     RandParamEsts2(:,j)=ParamEsts2';
%     RandGRho2(j)=tempModelRho2;
%     RandInitRho2(j)=tempModel2Rho2;
%     RandGLinR2(j)=tempR2FalsePosLinModel;
% end
% numRandDone=sum(~isnan(RandGNegLL));
% p_G=(sum(RandGNegLL(1:numRandDone)<=fvalMin)+1)/(numRandDone+1) %significance of bifurcation parameter modified by initial conditions (1-sided)
% p_IC=(sum(RandInitNegLL(1:numRandDone)<=fvalMin)+1)/(numRandDone+1) %significance of initial conditions (1-sided)
% p_Model_Lin=(sum(RandGRho2(1:numRandDone)-RandGLinR2(1:numRandDone)>=BestRho2-R2FalsePosLinModel)+1)/(numRandDone+1) %significance of model against linear regression (1-sided)
% p_d_G=(2*min(sum(RandParamEsts1(2,1:numRandDone)>=BestParamEsts(2)),sum(RandParamEsts1(2,1:numRandDone)<=BestParamEsts(2)))+1)/(numRandDone+1); %significance of parameter d aganist randomized bifurcation parameter g (2-sided)
% p_v_G=(2*min(sum(RandParamEsts1(4,1:numRandDone)>=BestParamEsts(4)),sum(RandParamEsts1(4,1:numRandDone)<=BestParamEsts(4)))+1)/(numRandDone+1); %significance of parameter v aganist randomized bifurcation parameter g (2-sided)
% p_d_IC=(2*min(sum(RandParamEsts2(2,1:numRandDone)>=BestParamEsts(2)),sum(RandParamEsts2(2,1:numRandDone)<=BestParamEsts(2)))+1)/(numRandDone+1); %significance of parameter d aganist randomized initial condition (2-sided)
% p_v_IC=(2*min(sum(RandParamEsts2(4,1:numRandDone)>=BestParamEsts(4)),sum(RandParamEsts2(4,1:numRandDone)<=BestParamEsts(4)))+1)/(numRandDone+1); %significance of parameter v aganist randomized initial condition (2-sided)
