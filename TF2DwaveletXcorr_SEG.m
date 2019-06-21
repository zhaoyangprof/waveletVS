clear all; clc; close all; warning off;

% add path
%addpath(genpath('/HRC_GPT_DATA1/zhaoyx/Software/matlab/SeismicLab/codes'));
addpath(genpath('C:\05_matlab\SeismicLab\codes'));
addpath scripts
flag = 'nnyn';

% load
% [csg,vsh] = readsegy('/HRC_GPT_DATA1/zhaoyx/Redatum/RedatumEchos/CRG/NoSurfaceConsisProc/ep1csg124.su');
%load('/HRC_GPT_DATA1/zhaoyx/Redatum/waveletRedatum/matlab/input/data/ep1trap65.mat');
[csg,vsh] = readsegy('ep1trap65.su');
shot = csg(1:1000,:);
tic;

wigb 
%
fig = 1;

%% define analysis parameters

wlen = 64;                        % window length (recomended to be power of 2)
dt = 0.002;
dj = 0.05;
s0 = 0.005;
%s0 = 0.01;
ub = 100;
lb = 1;
Nr = 7;
Nc = 3;
nt = 1000;
nti = 1000;
nfft = nti;
win = 25;
rec = 80;
FilterT = 12;
FilterF = 15;

t = linspace(0,dt*nt,nt);

%% create a 3D T-F-K mask for time-vary FK filtering

maskx = [1 40 80];
% masky = [11 76 11];
masky = [8 76 8];

[X,Y] = meshgrid(1:80,1:76);

mask = inpolygon(X,Y,maskx,masky);
mask(1:11,:) = 1;
maskeCube = repmat(mask,[1,1,nt]);
maskeCube(:,:,1:150) = 1;
maskeCubeFrotran = permute(maskeCube,[3,1,2]);
maskeCubeFrotranShift = fftshift(maskeCubeFrotran,3);
fileID1 = fopen('mask_1000_76_80_2000v.bin', 'w+');
fwrite(fileID1, maskeCubeFrotranShift, 'single');

maskeCube1 = permute(maskeCube,[1,3,2]);
maskeCubeShift = fftshift(maskeCube1,3);
fileID2 = fopen('mask_76_1000_80_2000v.bin', 'w+');
fwrite(fileID2, maskeCubeShift, 'single');

%% to generate a gaussian 40hz taget spetrum
x = [0 20 40 80 87 100 120 150 200 250]; x=x*2;
y = zeros(size(x)); y(2) = 6; y(3) = 8;
y1 = interp1(x,y,1:500,'spline');
ASR = abs([y1 fliplr(y1)]);
ASRFlag = 0;

down = zeros(nt,1);
down(1:250) = csg(1002:end,65);
down(75:end) = 0;

[wavedown, perioddown, scaledown, coidown, djdown,paramoutdown, kdown] = contwt(down,dt,[],dj,s0);

for iRec = 66
    
    up = csg(1:nt,iRec);
    
    conVS = xcorr(up,down);
    conVS = conVS(nt:end);
    
    
    if isequal(ASRFlag,1)
        % if doing ASR to downgoing wavefields
        ASRPha = angle(fft(down));
        % Replace 1D trace by ASR
        DaASR = ASR'.*complex(cos(ASRPha),sin(ASRPha));
        down = real(ifft(DaASR));
    end
    
    %% compute the CWT
    % [waveup, periodup, scaleup, coiup, djup,paramoutup, kup] = contwt(up,dt,[],dj,s0,[],'MORLET',6);
    [waveup, periodup, scaleup, coiup, djup,paramoutup, kup] = contwt(up,dt,[],dj,s0);
    
    %% Do the wavelet Cross-correlation
    wavex = zeros(length(periodup),2*nt-1);
    
    for ifreq = 1:length(perioddown)
        wavex(ifreq,:)  = xcorr(waveup(ifreq,:),wavedown(ifreq,:));
    end
    
    % do wavelet domain 2D filteringfi
    wavex2  = wavex(:,nti:end);
    amp = abs(wavex2);
    pha = angle(wavex2);
    
    % do expontional filtering in 2D amplitude domain
    % scale an scale back to original amplitude
    
    
    ampnew = exp(amp);
    
    %         for i = 1:nt
    %             maxtrc = max(abs(ampnew(:,i)))-1;
    %
    %             for j = 1:length(periodup)
    %                 if (ampnew(j,i)< maxtrc/3+1)
    %                     ampnew(j,i) = 1;
    %                 end
    %             end
    %
    %         end
    %         amp4 = log(ampnew);
    amp4 = amp;
    for it = FilterT+1:nt-FilterT
        %for it = 550:650
        for ifreq = FilterF+1:length(periodup)-FilterF
            WinW = ifreq-FilterF:ifreq+FilterF;
            WinL = it-FilterT:it+FilterT;
            maxtrc = max(abs(amp(WinW,WinL)));
            amp1 = amp(WinW,WinL)/maxtrc;
            amp2 = exp(amp1*4);
            amp2 = amp2-1;
            amp3 = amp2*maxtrc;
            amp4(WinW,WinL) = amp3;
            amp4(isnan(amp4)) = 0;
        end
    end
    
    % for it = FilterT+1:nt-FilterT
    %
    %     WinL = it-FilterT:it+FilterT;
    %     maxtrc = max(abs(amp(:,WinL)));
    %     amp1 = amp(:,WinL)/maxtrc*10;
    %     amp2 = exp(amp1);
    %     amp2 = amp2-1;
    %     amp3 = amp2*maxtrc;
    %     amp4(:,WinL) = amp3;
    %     amp4(isnan(amp4)) = 0;
    % end
    
    wavex1 = amp4.*complex(cos(pha),sin(pha));
    
    
    %% apply wavelet domain filter and sum over desired frequency
    %wavexrec = invcwt(wavex1, 'MORLET', scaleup, paramoutup, kup);
    
    
end


wavefgfilter = invcwt(wavex1, 'MORLET', scaleup, paramoutup, kup);


%%
wavefg = xcorr(up,down);
wavefg = wavefg(nt:end);

frequp = 1./periodup;

fig = fig + 1;
figure(fig);
plot(t,wavefg/max(abs(wavefg))*max(abs(wavefgfilter)));
hold on;
plot(t,wavefgfilter,'r');
xlim([1 2]);
ylim(1e-5*[-1.5 1.5]);
legend('cross-correlation','wavelet cross-correlation');
Name = '1D trace wavelet correlation comparison';
xlabel('Time (s)');
ylabel('Amplitude');
%title(Name);
saveFigure(fig,Name,flag);

%%
fig = fig + 1;
figure(fig);
subplot(4,1,[1:3]);
pcolor(t,frequp,abs(waveup));
shading flat;
ylabel('Frequency (Hz)');
ylim([5 100]);
set(gca,'Yscale','log');
set(gca,'YDir','reverse');
Name = 'Wavelet Domain Upgoing Wavefields';
title(Name);

subplot(4,1,4);
plot(t,up);
xlabel('Time (sec)');
ylabel('Amplitude');
saveFigure(fig,Name,flag);


fig = fig + 1;
figure(fig);
subplot(4,1,[1:3]);
pcolor(t,frequp,abs(wavedown));
shading flat;
ylabel('Frequency (Hz)');
set(gca,'YDir','reverse');
ylim([5 100]);
set(gca,'Yscale','log');
Name = 'Wavelet Domain Downgoing Wavefields';
title(Name);

subplot(4,1,4);
plot(t,down);
xlabel('Time (sec)');
ylabel('Amplitude');
saveFigure(fig,Name,flag);

close all;

%%
fig = fig + 1;
figure(fig);
subplot(4,1,[1:3]);
amp = amp/max(max(abs(amp)));
pcolor(t,frequp,amp);
shading flat;
h = colorbar;
ylabel(h, 'Magnitude')
ylabel('Frequency (Hz)');
%set(gca,'YDir','reverse');
ylim([5 100]);
%set(gca,'Yscale','log');
Name = 'Wavelet Cross-correlation before 2D TF filtering';
%title(Name);

subplot(4,1,4);
plot(t,wavefg);
xlabel('Time (sec)');
ylabel('Amplitude');
saveFigure(fig,Name,flag);

%%
fig = fig + 1;
figure(fig);
subplot(4,1,[1:3]);
amp4 = amp4/max(max(abs(amp4)));
pcolor(t,frequp,amp4);
h = colorbar;
ylabel(h, 'Magnitude')
shading flat;
ylabel('Frequency (Hz)');
%set(gca,'YDir','reverse');
ylim([5 100]);
%set(gca,'Yscale','log');
Name = 'Wavelet Cross-correlation after 2D TF filtering';
%title(Name);

subplot(4,1,4);
plot(t,wavefgfilter);
xlabel('Time (sec)');
ylabel('Amplitude');
%sublabel
saveFigure(fig,Name,flag);