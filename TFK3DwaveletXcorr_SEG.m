clear all; clc; close all; warning off;

% add path
%addpath(genpath('/HRC_GPT_DATA1/zhaoyx/Software/matlab/SeismicLab/codes'));
addpath(genpath('C:\05_matlab\SeismicLab\codes'));
addpath scripts
flag = 'nnnny';
fig = 1;
% load
% [csg,vsh] = readsegy('/HRC_GPT_DATA1/zhaoyx/Redatum/RedatumEchos/CRG/NoSurfaceConsisProc/ep1csg124.su');
%load('/HRC_GPT_DATA1/zhaoyx/Redatum/waveletRedatum/matlab/input/data/ep1trap65.mat');
[csg,vsh] = readsegy('ep1trap65.su');


%csg = squeeze(Cube(:,:,65));
shot = csg(1:1000,:);
tic;
%% define analysis parameters

wlen = 64;                        % window length (recomended to be power of 2)
dt = 0.002;
dj = 0.1;
%s0 = 0.008;
s0 = 0.0118;
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

%% create a 3D T-F-K mask for time-vary FK filtering

maskx = [1 40 80];
% masky = [11 76 11];
masky = [8 76 8];

[X,Y] = meshgrid(1:80,1:76);

mask = inpolygon(X,Y,maskx,masky);
mask(1:11,:) = 1;
%mask = flip(mask);
maskeCube = repmat(mask,[1,1,nt]);
%maskeCube(:,:,1:150) = 1;
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

for iRec = 1:rec
    
    up = csg(1:nt,iRec);
    
    tmp = xcorr(up,down);
    conVS(:,iRec) = tmp(nt:end);
    
    
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
    
    
    %     ampnew = exp(amp);
    %
    %     for i = 1:nt
    %         maxtrc = max(abs(ampnew(:,i)))-1;
    %
    %         for j = 1:length(periodup)
    %             if (ampnew(j,i)< maxtrc/2+1)
    %                 ampnew(j,i) = 1;
    %             end
    %         end
    %
    %     end
    %     amp4 = log(ampnew);
    
    % for it = FilterT+1:nt-FilterT
    %     %for it = 550:650
    %     for ifreq = FilterF+1:length(periodup)-FilterF
    %         WinW = ifreq-FilterF:ifreq+FilterF;
    %         WinL = it-FilterT:it+FilterT;
    %         maxtrc = max(abs(amp(WinW,WinL)));
    %         amp1 = amp(WinW,WinL)/maxtrc;
    %         amp2 = exp(amp1*2);
    %         amp2 = amp2-1;
    %         amp3 = amp2*maxtrc;
    %         amp4(WinW,WinL) = amp3;
    %         amp4(isnan(amp4)) = 0;
    %     end
    % end
    
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
    
    amp4 = amp;
    wavex1 = amp4.*complex(cos(pha),sin(pha));
    
    wavexCube(:,:,iRec) = wavex1;
    waveupCube(:,:,iRec) = waveup;
    
    %% apply wavelet domain filter and sum over desired frequency
    %wavexrec = invcwt(wavex1, 'MORLET', scaleup, paramoutup, kup);
    
    
end

wavedownCube = zeros(size(waveupCube));

for i = 1:10
    down = zeros(nt,1);
    down(1:250) = csg(1002:end,60+i);
    down(75:end) = 0;
    
    [wavedown, perioddown, scaledown, coidown, djdown,paramoutdown, kdown] = contwt(down,dt,[],dj,s0);
    
    wavedownCube(:,:,60+i) = wavedown;
    
end
% fileID1 = fopen('mask_76_1001_80', 'r+');
% tmp = fread(fileID1, 76*1000*80, 'single');
% maskeCubeShift = reshape(tmp,76,nt,80);

% TFK = fftshift(fft(wavexCube,[],3),3);
% TFKFilter = maskeCube.*TFK;
% wavexCubeFilter = ifft(ifftshift(TFKFilter,3),[],3);

TFKup = fft(waveupCube,[],3);
TFKdown = fft(wavedownCube,[],3);

TFK = fft(wavexCube,[],3);
TFKFilter = maskeCubeShift.*TFK;
wavexCubeFilter  = ifft(TFKFilter,[],3);

for iRec = 1:rec
    wavex2 = squeeze(wavexCubeFilter(:,:,iRec));
    tmp1  = invcwt(wavex2, 'MORLET', scaleup, paramoutup, kup);
    shotFilter(:,iRec) = tmp1;
end

%% TFK plotting
t = linspace(0,length(up)*dt,length(up));
f = 1./periodup;
%f = flip(f);
k = -2400:60:2400;
k = k(1:end-1);
k = 1./k;




%%
fig = fig + 1;
figure(fig);
%TFKdown(:,250:end,:) = 0;
TFKdownplot = abs(fftshift(TFKdown,3));
H = vol3d('CData',log10(TFKdownplot),'XData',t,'YData',f,'ZData',k);
caxis([-4 0]-1)
view([-30,30])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
zlabel('Wavenumber (1/m)')
%set(gca,'ydir','reverse');
%set(gca,'Yscale','log');
shading flat;
axis tight;
ylim([5 80]);
figureName = 'Downgoing Cube Before TFK Filtering';
saveFigure(fig,figureName, flag)

%%
fig = fig + 1;
figure(fig);
TFKupplot = abs(fftshift(TFKup,3));
H = vol3d('CData',log10(TFKupplot),'XData',t,'YData',f,'ZData',k);
caxis([-4 0]-3)
view([-30,30])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
zlabel('Wavenumber (1/m)')
%set(gca,'ydir','reverse');
%set(gca,'Yscale','log');
shading flat;
axis tight;
ylim([5 80]);
figureName = 'Upgoing Cube Before TFK Filtering';
saveFigure(fig,figureName, flag)



%%
close all;
fig = fig + 1;
figure(fig);
%subplot(2,1,1);
TFKplot = abs(fftshift(TFK,3));
%TFKplot1 = permute(TFKplot,[2,3,1]);
H = vol3d('CData',log10(TFKplot),'XData',t,'YData',f,'ZData',k);
% H = vol3d('CData',log10(TFKplot));
%xlim([0 45]);
caxis([-4 0]-4)
view([-30,30])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
zlabel('Wavenumber (1/m)')
%set(gca,'xdir','reverse');
%set(gca,'Yscale','log');
shading flat;
axis tight;
ylim([5 80]);
figureName = 'VS Cube Before TFK Filtering';

saveFigure(fig,figureName, flag)
% 
fig = fig+1;
figure(fig);

% TFKplotFilter = TFKplot.*maskeCube1;
TFKplotFilter = abs(fftshift(TFKFilter,3));
%TFKplotFilter1 = permute(TFKplotFilter,[2,3,1]);
H = vol3d('CData',log10(TFKplotFilter),'XData',t,'YData',f,'ZData',k);

caxis([-4 0]-4)
view([-30,30])
xlabel('Time (s)')
ylabel(' Frequency (Hz)')
zlabel('Wavenumber (1/m)')
%set(gca,'xdir','reverse');
%set(gca,'Yscale','log');
shading flat;
axis tight;
ylim([5 80]);
figureName = 'VS Cube After TFK Filtering';

saveFigure(fig,figureName, flag)

%%
fig = fig + 1;
figure(fig);
subplot(1,2,1);  
imagesc(1:rec,t,conVS);  
maxval = max(max(abs(conVS)));
colormap(gray);
ylabel('Time (s)');
xlabel('Receiver');
caxis(0.1*[-maxval maxval]);
%title('Before VS TFK filtering');

subplot(1,2,2);
imagesc(1:rec,t,shotFilter);
maxval = max(max(abs(shotFilter)));
colormap(gray);
xlabel('Receiver');
caxis(0.1*[-maxval maxval]);
%title('After VS TFK filtering');
figureName = 'Shot Gather TFK Filtering';
saveFigure(fig,figureName, flag)