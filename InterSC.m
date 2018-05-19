function InterSC()
%%
% Regev, Honey, Simony & Hasson, 2013
% Selective and Invariant Neural Responses to Spoken and Written Narratives
%
% Shang, Linlin
%
% State Key Laboratory of Cognitive Neuroscience and Learning-BNU
% Xinjiekouwai Street 19, Beijing Haidian 100875
%
% Email: shangll@outlook.com | shangll@bnu.edu.cn
%
% URL: brdl.bnu.edu.cn/

%%
clear all;

% input data
% enter work directory
groupA = inDat('D:/Usr/Documents/BNU/A');
groupB = inDat('D:/Usr/Documents/BNU/B');

% raw corr coeffiencent
% both are 3D
wRA = withinCorr(groupA);
wRB = withinCorr(groupB);
aRAB = acrossCorr(groupA, groupB);

thresholdA = getWithinThreshold(groupA);
thresholdB = getWithinThreshold(groupB);

thresholdAB = getAcrossThreshold(groupA, groupB);

% activation
wRA(wRA<=thresholdA) = 0;
wRA(wRA>thresholdA) = 1;
save('wRA.mat', 'wRA');

wRB(wRB<=thresholdB) = 0;
wRB(wRB>thresholdB) = 1;
save('wRB.mat', 'wRB');

aRAB(aRAB<=thresholdAB) = 0;
aRAB(aRAB>thresholdAB) = 1;
save('aRAB.mat', 'aRAB');

end

%% imput & classify data

function group = inDat(filePath)
files = dir(fullfile(filePath,'*.mat'));

for n = 1:length(files)
    targetFile = fullfile(filePath, files(n).name);
    dat = cell2mat(struct2cell(load(targetFile)));
    group(:,:,:,:,n) = dat;
end

end
%% correlation within a condition

% groupX is a 5D matrix
function withinRX = withinCorr(groupX)
[ro, col, pg, t, chk] = size(groupX);

% 4D time series
allSubj = sum(groupX, 5);

for subjN = 1:chk
    residues(:,:,:,:,subjN) = (allSubj - groupX(:,:,:,:,subjN))/(chk-1);
    
    for z = 1:pg
        for x = 1: ro
            for y = 1:col
                r = corrcoef(squeeze(groupX(x,y,z,:,subjN)), squeeze(residues(x,y,z,:,subjN)));
                allR(x, y, z, subjN) = r(1,2);
                
            end
        end
    end
    
    % 4D to 3D matrix
    withinRX = mean(allR, 4);
 
end

end


%% correlation of across conditions

function acrossR = acrossCorr(groupX, groupY)
[Xro, Xcol, Xpg, Xt, Xchk] = size(groupX);

% the average time series of all the subj in the other group
% 4D
meanY = mean(groupY, 5);

for subjX = 1:Xchk
    for z = 1:Xpg
        for x = 1:Xro
            for y = 1:Xcol
                r = corrcoef(squeeze(groupX(x,y,z,:,subjX)), squeeze(meanY(x,y,z,:)));
                allX(x, y, z, subjX) = r(1,2);
                
            end
        end
    end
   
end

acrossR = mean(allX, 4);

end

%% Fourier Transf
% References:
% Theiler, Galdrikian, Longtin, Eubank & Farmer, 1992
% Theiler, Eubank, Galdrikian, Longtin & Farmer, 1992
%
% Alexandros Leontitsis
% me00743@cc.uoi.gr | leoaleq@yahoo.com
function newDat = fourierRandInv(dat)
[ro, col, pg, t, chk] = size(dat);
for subj = 1:chk
    for z = 1:pg
        for x = 1:ro
            for y = 1:col
                % FFT
                rawDat = squeeze(dat(x,y,z,:,subj));
                Fgr=fft(rawDat);
                m=abs(Fgr);
                p=angle(Fgr);
                i=sqrt(-1);
                h=floor(t/2);
                if rem(t,2)==0
                    p1=rand(h-1,1)*2*pi;
                    p(2:t)=[p1' p(h+1) -flipud(p1)'];
                    m=[m(1:h+1);flipud(m(2:h))];
                else
                    p1=rand(h,1)*2*pi;
                    p(2:t)=[p1 -flipud(p1)];
                end
                
                % 5D
                newDat(x,y,z,:,subj)=m.*exp(i*p);
                newDat(x,y,z,:,subj)=real(ifft(newDat(x,y,z,:,subj)));
            
            end
        end
     end
end

end

%% get within threshold

function WR = getWithinThreshold(group)
for repN = 1:1000
    % fft
    newGroup = fourierRandInv(group);
    % group & newGroup ¡ú 5D
    
    % set up null distribution
    nullDistr(:,:,:,repN) = withinCorr(newGroup);

end

% find maximum values for each voxel
[d1, d2, d3, repN] = size(nullDistr);
for z = 1:d3
    for x = 1:d1
        for y = 1:d2
            repR = squeeze(nullDistr(x,y,z,:));
            maxR(x,y,z) = max(repR);
            
        end
    end
end

% find the value at the top 5% of the distribution
tempV = sort(maxR(:));
n = floor(.95*length(tempV));
WR = tempV(n);

end

%% get across threshold

function AR = getAcrossThreshold(groupX, groupY)
for repN = 1:1000
    % fft
    newGroupX = fourierRandInv(groupX);
    newGroupY = fourierRandInv(groupY);
    % group & newGroup ¡ú 5D
    
    % set up null distribution
    nullDistr(:,:,:,repN) = acrossCorr(newGroupX, newGroupY);

end

% find maximum values for each voxel
[d1, d2, d3, repN] = size(nullDistr);
for z = 1:d3
    for x = 1:d1
        for y = 1:d2
            repR = squeeze(nullDistr(x,y,z,:));
            maxR(x,y,z) = max(repR);
            
        end
    end
end

% find the value at the top 5% of the distribution
tempV = sort(maxR(:));
n = floor(.95*length(tempV));
AR = tempV(n);

end