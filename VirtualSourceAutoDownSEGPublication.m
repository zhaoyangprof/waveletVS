% perform virtual source for synthetic
% =========================================================================
%
clear; clc; close all; warning('off','all');


addpath(genpath('C:\05_matlab\SeismicLab\codes'));
addpath scripts
flag = 'nnyn';


datapath = 'input/data/';
hdrpath ='input/hdr/';


% load conventional redatum data
load input/VSCubeCorr.mat
RedatumCube = VSCube;
clear VSCube;
% clear up the old data
! rm output/ASRVS1/*
! rm *.eps *.png *~

% to generate a gaussian 40hz taget spetrum
x = [0 20 40 80 87 100 120 150 200 250]; x=x*2;
y = zeros(size(x)); y(2) = 6; y(3) = 8;

y1 = interp1(x,y,1:500,'spline');
ASR = abs([y1 fliplr(y1)]);

% set up paramters/
dt = 0.002;
tcut = 596/2;
nti = 1000 ;
rec = 80;
ep = 1;
aperture = 2400;
mute =60;
fig = 1;
f0 = 40;
depth = 30;
Inline = 9;
t=0:dt:nti*dt;
w = hamming((aperture-1)*2);
taperold = [ones(1,1);w(end/2+1:end)];
tapernew = taperold.^20;

v = hamming((207*2));
timetaper = v(end/2+1:end);

% set up a perfect ricker wavelet center at 40Hz
[rickerwave,temp] = ricker(f0,dt);
ricker = zeros(250,1);
ricker(1:length(rickerwave)) = rickerwave;

% read in data
showtime = linspace(0,2,nti);
showfreq = linspace(0,500,nti);
scsg = zeros(nti,rec,rec);



% define target window
tarwin = 450:950;
tarwinLength = length(tarwin);

pickedvs = 60;
%% Assign ASR and FK ASR into different Da and save for cross-correlation VS
tic
%for iEp = 1:5:13
for iEp = [1 11]
    tic
    iEp
    
    dataname = strcat('input/vs',num2str(iEp),'crg_data.mat')
    hdrname  = strcat('input/vs',num2str(iEp),'crg_header.mat')
    
    load(dataname)
    load(hdrname)
    
    nt = size(D,1);
    nshot = size(D,2)/rec;
    
    Cube = reshape(D,nt,rec,nshot);
    HDR = reshape(H,rec,nshot);
    
    
    %% perform virtual source redatuming
    % loop over all shotlines and each shotline order
    % VSCube = zeros(nti,rec,Shots);
    VSCubeNew = zeros(nti,rec,rec);
    VSCubeOld = zeros(nti,rec,rec);
    
    % Everthing based on each vritual source as a unit
    %for ivs = 20:20:60
    for ivs = pickedvs;
        
        ivs
        %% Form Da cube
        % Only extract ivs number of Da from Da Cube. other Da are trash
        DaCRG = squeeze(Cube(1001:end,:,:));
        DaCRG(201:end,:,:) = 0;
        %DaFK  = fft(DaCRG,nti*2);
        
        %% Reshape to be 3D ASR Cube by Inline shotlines order
        % load Rl Cube
        RlCRG = Cube(1:nti,:,:);
        %RlFK  = fft(RlCRG,nti*2);
        
        %% fft the initial redatum cube
        %RedatumFK = fft(RedatumCube,nti*2);
        
        cnt = 1;
        
        for iShot = 1:nshot
            
            % Apply spatial taper in time domain
            x = (HDR(ivs,iShot).sx/10-HDR(ivs,iShot).gx/10)^2;
            y = (HDR(ivs,iShot).sy/10-HDR(ivs,iShot).gy/10)^2;
            dist = sqrt(x+y);
            
            if dist < length(tapernew)-1
                
                offsetaper = tapernew(round(dist+1));
                offsetaperold = taperold(round(dist+1));
                
                downgoing = DaCRG(:,ivs,iShot);
                
                
                
                [downgoingbp] =  bp_filter(downgoing,0.002,12,15,65,70);
                downfft = abs(fft(xcorr(downgoingbp)))/max(abs(fft(xcorr(downgoingbp))));
                downcoeff = sum(downfft);%/length(downfft);
                
                % band-pass downgoing
                
                downgoing1 = downgoingbp*downcoeff*offsetaperold;%*offsetaper;
                downgoing2 = downgoing*offsetaper;
                
                % output for SEG publication
                if isequal(ivs,ivs)
                    
                    if iEp < 7
                        Xcube(cnt) = HDR(ivs,iShot).sx/10;
                        Ycube(cnt) = HDR(ivs,iShot).sy/10;
                    else
                        Xcube(cnt) = HDR(ivs,iShot).sx;
                        Ycube(cnt) = HDR(ivs,iShot).sy;
                    end
                    offsetaperCube(cnt) = offsetaper;
                    downcoeffCube(cnt) = downcoeff;
                    
                    cnt = cnt + 1;
                    
                end
                
                tmp3 = zeros(nti,rec);
                tmp4 = zeros(nti,rec);
                
                for iRec = 1:rec
                    
                    
                    % ! Check the source location for direct and reflection data, and then xcorr
                    % if ((sx_da==sx_rl(ig)) .and. (sy_da==sy_rl(ig))) then
                    
                    if isequal(HDR(ivs,iShot).sx,HDR(iRec,iShot).sx) && isequal(HDR(ivs,iShot).sy,HDR(iRec,iShot).sy)
                        
                        
                        redatum = RedatumCube(:,iRec,ivs);
                        upgoing = RlCRG(:,iRec,iShot);
                        
                        % perform conv back to upgoing
                        if dist < 5
                            redtmp = conv(downgoing,redatum);
                            redupgather(:,iRec) = redtmp(1:nti);
                        end
                        
                        
                        %uptarget = upgoing(tarwin);
                        %reduptarget = redupgoing(tarwin);
                        
                        % get the matched filter
                        %[mfilt,tm]=matchf(reduptarget,uptarget,t,50,2);
                        
                        
                        % perform one side (primary) redatuming;;
                        tmp1 = xcorr(downgoing1,upgoing);
                        tmp2 = xcorr(downgoing2,upgoing);
                        
                        % apply matched filter
                        %tmp1 = conv(tmp(1:nti),mfilt);
                        %tmp2=tmp1(251:end-250);
                        tmp3(:,iRec) = tmp1(1:nti);
                        tmp4(:,iRec) = tmp2(1:nti);
                        
                        %sum(coffweight);
                        % taper(round(dist));
                        
                        % clear out NaN value
                        if isnan(tmp3(:,iRec))
                            tmp3(:,iRec) = 0.0;
                            tmp4(:,iRec) = 0.0;
                        end
                        
                        
                        
                        
                    end
                    
                end
                
                tmp31 = flipud(tmp3(1:nti,:));
                tmp41 = flipud(tmp4(1:nti,:));
                
                
                VSCubeNew(:,:,ivs) = tmp31 + VSCubeNew(:,:,ivs);
                VSCubeOld(:,:,ivs) = tmp41 + VSCubeOld(:,:,ivs);
                
                if dist < 5
                    upgather = RlCRG(:,:,iShot);
                end
                
                
            end
            
            
            
        end
        
        WeightName = strcat('WeightTaper_vs',num2str(ivs),'ep',num2str(iEp));
        save(WeightName,'Xcube','Ycube','offsetaperCube','downcoeffCube');
        
        %% save all VS files
        %         binname = strcat('output/ASRVS1/','ep',num2str(iEp),'vs',num2str(ivs),'.bin');
        %         suname = strcat('output/ASRVS1/','ep',num2str(iEp),'vs',num2str(ivs),'.su');
        %         fileID = fopen(binname,'w+');
        %         fwrite(fileID, VSCubeNew(:,:,ivs), 'single');
        %         fclose(fileID);
        %
        %         binname = strcat('output/ASRVS1/','ep',num2str(iEp),'vs',num2str(ivs),'_old.bin');
        %         suname = strcat('output/ASRVS1/','ep',num2str(iEp),'vs',num2str(ivs),'_old.su');
        %         fileID = fopen(binname,'w+');
        %         fwrite(fileID, VSCubeOld(:,:,ivs), 'single');
        %         fclose(fileID);
        
        %% save into SU files
        %         setenv('bin',binname);
        %         setenv('su',suname);
        %         !suaddhead ns=1000 <$bin | sushw key=dt a=2000 >$su
        
        %% save for taper information
        close all;
        %% plot adaptive weightin
        %
        fig = fig + 1;
        
        load WeightTaper_vs60ep1.mat
        ivs = 60;
        
        Xcube1 = Xcube - min(Xcube);
        Ycube1 = Ycube - min(Ycube);
        
        if iEp < 7
            vsx1 = HDR(ivs,iShot).gx/10;
            vsy1 = HDR(ivs,iShot).gy/10;
        else
            vsx1 = HDR(ivs,iShot).gx;
            vsy1 = HDR(ivs,iShot).gy;
        end
        
        vsx = vsx1 - min(Xcube);
        vsy = vsy1 - min(Ycube);
        
        
        offsetaperCube1 = offsetaperCube-min(offsetaperCube);
        offsetaperCube0 = offsetaperCube1.^0.9;
        offsetaperCube2 = offsetaperCube0/max(offsetaperCube0);
        
        
        figure(fig);
        
        subplot(3,1,1);
        scatter(Xcube1,Ycube1,2,offsetaperCube2,'fill');
        hold on;
        plot(vsx,vsy,'v','MarkerFaceColor','k');
        hold off;
        h = colorbar;
        ylabel(h, 'Magnitude')
        grid on;
        xlabel('X (m)'); ylabel('Y (m)');
        %axis square;
        name = strcat('convention_vs_autodown_taper_vs',num2str(ivs),'ep',num2str(iEp));
        
        
        
        %downcoeffCube1 = downcoeffCube.*offsetaperCube2;
        downcoeffCube1 = downcoeffCube.^0.7;
        downcoeffCube2 = downcoeffCube1/max(downcoeffCube1);
        subplot(3,1,2)
        scatter(Xcube1,Ycube1,2,downcoeffCube2,'fill');
        hold on;
        plot(vsx,vsy,'v','MarkerFaceColor','k');
        hold off;
        h = colorbar;
        ylabel(h, 'Magnitude')
        grid on;
        xlabel('X (m)'); ylabel('Y (m)');
        %axis square;
        
        load WeightTaper_vs60ep11.mat
        Xcube1 = Xcube - min(Xcube);
        Ycube1 = Ycube - min(Ycube);
        
        downcoeffCube3 = downcoeffCube.^0.7;
        downcoeffCube4 = downcoeffCube3/max(downcoeffCube3);
        subplot(3,1,3)
        scatter(Xcube1,Ycube1,2,downcoeffCube4,'fill');
        hold on;
        plot(vsx,vsy,'v','MarkerFaceColor','k');
        hold off;
        h = colorbar;
        ylabel(h, 'Magnitude')
        grid on;
        xlabel('X (m)'); ylabel('Y (m)');
        %axis square;
%         
%         saveFigure(fig, name, flag);
%         
%         %clear downcoeffCube offsetaperCube
    end
    toc
    
    
    
end
toc


