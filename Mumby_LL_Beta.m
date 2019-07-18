function LL = Mumby_LL_Beta(g,avgAllPt,firstPt,firstM,tEnd,params)
digits(5); %decrease vpa precision

%g, grazing rate on macroalgae and algal turf (array of values)

r=params(1); %coral overgrowth on algal turf
d=params(2); %coral mortality
%a=params(3); %macroalgae overgrowth rate on coral
v=params(4); %gamma, macroalgae overgrowth rate on algal turf
a=v; %macroalgae overgrowth rate on coral - set to be identical to gamma (v)
phi=params(5); %precision parameter for coral cover (large when dispersion is relatively low) = omega+tau (beta's shape parameters); >2

%known solutions for coral cover
%known solutions for coral cover
Sol1=0*g;
% Sol2=(r-d+(r+v)*(d*v-r*v-v^2+((v*(v*d^2+2*d*r*v+2*d*v^2+v*r^2-4*g*r^2+2*r*v^2-4*g*r*v+v^3)).^(.5)))/(2*(v^2+r*v)))/r;
% Sol2(imag(Sol2)~=0)=-1;
% Sol2_M=1/2+((d*(d+2*r+2*v)+r*(r-4*g*r/v+2*v-4*g)+v^2).^.5-d)/(2*(r+v));
% Sol2_M(imag(Sol2_M)~=0)=NaN;
Sol2_C_M=[[-(d - r + ((a + r)*(a*r^2 - 2*a^2*d - r*(a^4 + 2*a^3*r - 4*g*a^3 + 2*a^2*d*v + a^2*r^2 - 8*g*a^2*r + 4*g*a^2*v + 2*a*d*r*v - 4*g*a*r^2 + 4*g*a*r*v + d^2*v^2).^(1/2) + a^2*r - 2*a*d*r + 2*a*d*v + d*r*v))/(2*(a^3 + 2*a^2*r - v*a^2 + a*r^2 - v*a*r)))/r]; [(a*r^2 - 2*a^2*d - r*(a^4 + 2*a^3*r - 4*g*a^3 + 2*a^2*d*v + a^2*r^2 - 8*g*a^2*r + 4*g*a^2*v + 2*a*d*r*v - 4*g*a*r^2 + 4*g*a*r*v + d^2*v^2).^(1/2) + a^2*r - 2*a*d*r + 2*a*d*v + d*r*v)/(2*(a^3 + 2*a^2*r - v*a^2 + a*r^2 - v*a*r))]];
Sol2=Sol2_C_M(1,:);
Sol2(imag(Sol2)~=0)=-1;
Sol2_M=Sol2_C_M(2,:);
Sol2_M(imag(Sol2_M)~=0)=NaN;
Sol3=0*g+1-d/r;
% Sols=[Sol1;real(Sol2);Sol3];
%
% for i=1:length(g)
%     if Sol2(i)>Sol3(i)
%         Sols(3,i)=NaN;
%         Sols(2,i)=NaN;
%     end
%     if Sol2(i)<Sol1(i)
%         Sols(1,i)=NaN;
%         Sols(2,i)=NaN;
%     end
% end

% %restrict avgAllPt to min 0.01 and max 0.99:
avgAllPt(avgAllPt<=0)=0.001;
avgAllPt(avgAllPt>=1)=0.999;


%evaluate fit of data to modelled equilibria
LLmeans=[]; %store log likelihoods
X=[];
for i=1:length(avgAllPt)
    %parameterized dynamic system (for figuring out basin of attraction)
    fMumby = @(t,X) [X(1)*(r*(1-X(1)-X(2))-d-a*X(2)); X(2)*(a*X(1)-g(i)/(X(2)+(1-X(1)-X(2)))+v*(1-X(1)-X(2)))]; %X(1) is coral cover, X(2) is macroalgal cover
    [~,X] = ode23(fMumby,[0 tEnd],[firstPt(i) firstM(i)]);     % Runge-Kutta 2nd/3rd order ODE solver
    
    [omega1,tau1]=Beta_Params(Sol1(i),phi);
    [omega2,tau2]=Beta_Params(Sol2(i),phi);
    [omega3,tau3]=Beta_Params(Sol3(i),phi);
    if Sol2(i)>Sol3(i) %solution is 0
        LLmeans(i)=log(betapdf(avgAllPt(i),omega1,tau1));
    elseif Sol2(i)<Sol1(i) %solution is the higher equilibrium
        LLmeans(i)=log(betapdf(avgAllPt(i),omega3,tau3));
    else %solution depends on initial condition
        % if firstPt(i)/firstM(i)<Sol2(i)/Sol2_M(i) %initial conditions below the linearly approximated separatrix, thus converges on 0 
        if abs(X(end,1)-Sol1(i))<abs(X(end,1)-Sol3(i))  %bistable region, projected long-run state close to lower stable state, below unstable solution
            LLmeans(i)=log(betapdf(avgAllPt(i),omega1,tau1));
            % elseif firstPt(i)/firstM(i)==Sol2(i)/Sol2_M(i) %on unstable point, stays there
            %    LLmeans(i)=log(betapdf(avgAllPt(i),omega2,tau2));
        else %coral cover converges on the higher stable equilibrium
            LLmeans(i)=log(betapdf(avgAllPt(i),omega3,tau3));
        end
    end
    X=[];
end
clearvars X
LL=sum(-LLmeans); %sum of -log likelihoods