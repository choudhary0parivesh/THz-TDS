clear all;
close all;
clc;

samname = 'Sample.mat'; %Sample File
refname = 'Reference.mat'; %Reference File 

%Loading the reference and collecting retval data (already avaraged)
load(refname);
amprefx = ydataX;
amprefy = ydataY;
load(samname);
ampx = ydataX_s; %X-voltage value from lock-in
ampy = ydataY_s; %Y-voltage value from lock-in
t = xdataX;% Translation stage position converted time


ampmaxstor = -10000;
for theta =-pi:.001:pi
    ampnew = cos(theta).*ampx + sin(theta).*ampy;
    [ampmax,ix] = max(ampnew);
    if (ampmax > ampmaxstor) 
    ampmaxstor = ampmax;
    imax = ix;
    ampbest = ampnew;
    thetamax = theta;
    end;
end;
amp = ampbest;
% amp = fliplr(amp')';
thetamaxsample=thetamax;

ampmaxstor = -10000;
for theta =-pi:.001:pi
    ampnew = cos(theta).*amprefx + sin(theta).*amprefy;
    [ampmax,ix] = max(ampnew);
    if (ampmax > ampmaxstor) 
    ampmaxstor = ampmax;
    imax = ix;
    ampbest = ampnew;
    thetamax = theta;
    end;
end;
ampref = ampbest;
% amp3 = fliplr(amp3')';
thetamaxsample=thetamax;

% z = ; %Translation stage position



[refxx,refyy] = max(ampref);
deltaref = refyy;
[samxx,samyy] = max(amp);
deltasam = samyy;
delta_abs = (refyy-samyy)*(t(1)-t(2));

[Max,mref]=max(ampref); % Placing the peak of the reference measurement as the zero point in time domain.
tmax=t(mref);
t = t - tmax;

% Create the time domain signal figure
figure
plot(t,ampref,'r','linewidth',1.5)
hold on
plot(t,amp,'k','linewidth',1.5)
xlabel('Time [ps]')
ylabel('Electric field strength [V/m]')
legend('Air','Silicon')


%Setting parameters
N=length(t);                          % Number of points
tstep = t(2)-t(1);                    % Time-step      
fstep = 1./(N.*tstep);            % Frequency-step
freq = fstep.*((.5-N/2):((N-.5)/2));


%%% - FFT %%%
samfft = fftshift(fft(amp,N)) ;
reffft = fftshift(fft(ampref,N));

figure
semilogy(freq ,abs(reffft),'r','linewidth',1.5)
hold on
semilogy(freq ,abs(samfft),'k','linewidth',1.5)
legend('Air','Si')
xlim([0 max(freq)])
% ylim([10^-5 1.05])
xlabel('Frequency [THz]')
ylabel('Spectral Amplitude [a.u.]')


freq=freq(length(samfft)/2:length(samfft));
Hx_i = samfft(length(samfft)/2:end)./reffft(length(samfft)/2:end);
Hx=abs(Hx_i);

%%%%%%%%%%%%%%%%%%%%%%%%

%Calculating n and K

%%%%%%%%%%%%%%%%%%%%%%%%%
% DELTA_T =  (2*(ref(1,3)-sample(1,3))./.3)
% Hx_THETA = unwrap(angle(Hx_i))-2*pi*DELTA_T*freq';

Hx_THETA = unwrap(angle(Hx_i))-2*pi;

figure
plot(freq,Hx_THETA,'r','linewidth',1.5)
xlim([0 10])
xlabel('Freq - THz')

n_air = 1.00027;
thick = 0.665*1e-3;% thickness of sample in meters
n_sample = zeros(1,length(Hx_i));
c = 3e8;

for j= 1:length(Hx_THETA)
   n_sample(j) =  n_air - ((Hx_THETA(j))*c)./(2*pi*freq(j)*thick*1e12);
   %kappa_sample(j) = -1e-12*(c/(2*pi*freq(j)*thick))*log(Hx(j)*((n_air+n_sample(j))^2/(4*n_air*n_sample(j)))); 
   alfa_sample(j) = -log(Hx(j)/(1-(n_sample(j)-1)^2/(n_sample(j)+1)^2))/(100*thick);
end



figure
plot(freq,n_sample,'r','linewidth',1.5)
hold on
plot(freq,ones(length(n_sample))*3.41,'k')
xlim ([0 2.5])
ylim ([3.1 3.5]) 
xlabel('Frequency [THz]')
ylabel('refractive index')



figure
plot(freq,alfa_sample,'r','linewidth',1.5)
xlim ([0 3])
% ylim ([35 70])
xlabel('Frequency [THz]')
ylabel('absorption Coefficent  (1/cm)')
x(1:length(reffft))=1;


% % ALFAMAX
% 
reffft = reffft./std(reffft(end/1.5:end));
alpha_max = 2*log(abs(reffft(length(reffft)/2:end)).*(4*n_sample./(n_sample+1).^2));

figure
plot(freq,(alpha_max)*10^-2/thick,'r','linewidth',1.5)
hold on
plot(freq,(alfa_sample),'b','linewidth',1.5)
legend('alfa_{max}','alfa_{measured}')
xlim ([0 max(freq)])
xlabel('Frequency [THz]')
ylabel('CM^{-1}')






