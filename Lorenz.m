%Lorenz.m
%Edward Tekwa
%Compute Lorenz system solutions

TimeData=string(datetime);

sigma = 10;
beta = 8/3;
rho = 28;
tEnd = 400;
initXYZ=[1 1 1];
tRec=ceil(tEnd*0.25)+1; % time after which to record random points
numPts = 10; %number of extrema to record and plot per time-series
%numAvgPts = 100; %number of random points to use for averaging at the end of time-series
XYZindex=3; %set to vary and plot initial X, Y, or Z (1, 2, or 3)
%quant=numPts/(tEnd-tRec+1);
labels=['X','Y','Z'];

rhoRange=[0:1:45];
XinitRange=[-60:.25:60];
%XinitRange=[0:5:60];
XYZs=zeros(length(XinitRange),length(rhoRange),numPts,3);
XYZAvgs=zeros(length(XinitRange),length(rhoRange),3);


for i=1:length(rhoRange)
    rho=rhoRange(i)
    for j=1:length(XinitRange)
        initXYZ(1)=XinitRange(j);
        f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
        [t,a] = ode45(f,[0 tEnd],initXYZ);     % Runge-Kutta 4th/5th order ODE solver
        [d,ix]=min(abs(t-[tRec:tEnd])); %find indices (ix) in time series corresponding to integer time steps
        %XYZs(j,i,:,:)=a(datasample(ix,numPts,'Replace',false),:); %should use near integer times instead of random sample time-steps
        XYZtemp=a(ix,:);
        diffXYZtemp=diff(XYZtemp);
        sortedDiffXYZtemp=sort(sum(abs(diffXYZtemp),2));
        candidateXYZs=XYZtemp(find(sum(abs(diffXYZtemp),2)<=sortedDiffXYZtemp(numPts+1))+1,:);
        XYZs(j,i,:,:)=candidateXYZs(1:numPts,:);
        %XYZAvgs(j,i,:)=mean(XYZs(j,i,:,:));
        XYZAvgs(j,i,:)=mean(XYZtemp);
    end
end


set(0,'defaultAxesFontSize',16)
scrsz = get(0,'ScreenSize');
fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3)/1 scrsz(4)/1]);
colormap jet
jitterSize=5;
for index=1:3
    XYZindex=index;
    subplot(2,3,index)
    %scatter(repmat(reshape(repmat(rhoRange,length(XinitRange),1),1,[]),1,numPts),reshape(XYZs(:,:,:,XYZindex),1,[]),0.1,reshape(repmat(XinitRange,length(rhoRange)*numPts,1)',1,[]),'.');
    scatter(repmat(reshape(repmat(rhoRange,length(XinitRange),1),1,[]),1,numPts),reshape(XYZs(:,:,:,XYZindex)+jitterSize*(rand(size(XYZs(:,:,:,XYZindex)))-0.5),1,[]),0.1,reshape(repmat(XinitRange,length(rhoRange)*numPts,1)',1,[]),'.');
    xlabel '\rho'
    ylabel (labels(index))
    %colorbar
    xlim([0 rhoRange(end)])
    
    subplot(2,3,index+3)
    scatter(reshape(repmat(rhoRange,length(XinitRange),1),1,[]),reshape(XYZAvgs(:,:,XYZindex),1,[]),0.1,reshape(repmat(XinitRange,length(rhoRange),1)',1,[]),'.');
    xlabel '\rho'
    ylabel (['mean(' labels(index) ')'])
    %colorbar
    xlim([0 rhoRange(end)])
end

fig=figure ('Color', [1 1 1],'Position',[1 scrsz(2) scrsz(3) scrsz(4)/3]);
colormap jet
rhoFocals=[5 10 15 20 23 24 24.5 28];
for i=1:8
    subplot(2,4,i)
    rhoFocal=find(rhoRange==rhoFocals(i));
    scatter(XinitRange,XYZAvgs(:,rhoFocal,1),24,XinitRange,'.');
    xlim([min(XinitRange) max(XinitRange)])
    ylim([-10 10])
    title (['\rho=' num2str(rhoFocals(i))])
end

FoodWebFile=sprintf('Lorenz_%s.mat', TimeData);
save(FoodWebFile);