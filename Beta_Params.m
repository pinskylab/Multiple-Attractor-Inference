function [omega,tau] = Beta_Params(mode,phi)

%Edward Tekwa Feb 25, 19
%inputs: mode (model stable state); phi (precision=omega+tau)
%outputs: omega, tau (beta's shape parameters)

%mode=(omega-1)/(tau+omega-2)=(omega-1)/(phi-2); omega=mode*(phi-2)+1
%***assumes there is one mode (no over-dispersion)
omega=mode*(phi-2)+1;
%phi=tau+omega;
%tau=phi-omega=phi-mode*(phi-2)-1
tau=phi-mode*(phi-2)-1;