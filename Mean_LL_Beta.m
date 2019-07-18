function LL = Mean_LL_Beta(avgAllPt,params) %fit beta mode and precision to data regardless of g
digits(5); %decrease vpa precision

%g, grazing rate on macroalgae and algal turf (array of values)

mode=params(1); %coral overgrowth on algal turf
phi=params(2); %coral mortality

% %restrict avgAllPt to min 0.01 and max 0.99:
avgAllPt(avgAllPt<=0)=0.01;
avgAllPt(avgAllPt>=1)=0.99;

%evaluate fit of data to mode
LLmodes=[]; %store log likelihoods
for i=1:length(avgAllPt)
    [omega,tau] = Beta_Params(mode,phi);
    LLmodes(i)=log(betapdf(avgAllPt(i),omega,tau));
end
LL=sum(-LLmodes); %sum of -log likelihoods