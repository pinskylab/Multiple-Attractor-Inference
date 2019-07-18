function SS = Mumby_SS(g,avgAllPt,firstPt,firstM,tEnd,params)
digits(5); %decrease vpa precision
%returns sum of squares

%g, grazing rate on macroalgae and algal turf (array of values)

r=params(1); %coral overgrowth on algal turf
d=params(2); %coral mortality
%a=params(3); %macroalgae overgrowth rate on coral
v=params(4); %gamma, macroalgae overgrowth rate on algal turf
a=v;

%known solutions for coral cover
Sol1=0*g;
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

%evaluate fit of fishery mean harvest rate data to modelled equilibria
SSmeans=[];
X=[];
for i=1:length(avgAllPt)
    %simulate to decide basin of attraction
    fMumby = @(t,X) [X(1)*(r*(1-X(1)-X(2))-d-a*X(2)); X(2)*(a*X(1)-g(i)/(X(2)+(1-X(1)-X(2)))+v*(1-X(1)-X(2)))]; %X(1) is coral cover, X(2) is macroalgal cover
    [~,X] = ode23(fMumby,[0 tEnd],[firstPt(i) firstM(i)]);     % Runge-Kutta 2nd/3rd order ODE solver
    
    if Sol2(i)>Sol3(i) %solution is 0
        SSmeans(i)=(avgAllPt(i)-Sol1(i))^2;
    elseif Sol2(i)<Sol1(i) %solution is the higher equilibrium
        SSmeans(i)=(avgAllPt(i)-Sol3(i))^2;
    else %solution depends on initial condition
%        if firstPt(i)/firstM(i)<Sol2(i)/Sol2_M(i) %first known coral cover is less than the unstable equilibrium, thus converges on 0
        if abs(X(end,1)-Sol1(i))<abs(X(end,1)-Sol3(i))
            SSmeans(i)=(avgAllPt(i)-Sol1(i))^2;
%         elseif firstPt(i)/firstM(i)==Sol2(i)/Sol2_M(i) %on unstable point, stays there
%             SSmeans(i)=(avgAllPt(i)-Sol2(i))^2;
        else %coral cover converges on the higher stable equilibrium
            SSmeans(i)=(avgAllPt(i)-Sol3(i))^2;
        end
    end
    X=[];
end
clearvars X
SS=sum(SSmeans);