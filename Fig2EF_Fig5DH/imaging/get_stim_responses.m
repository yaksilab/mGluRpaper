function responses = get_stim_responses(data,preBuff,postBuff,pix2micron)

% get stim indices 
assert(size(unique(vertcat(data.calStims), 'rows'),1) == 1) % ensure stim indices are the same 
stimIdx = data(1).calStims;

% go through each recording
patchedResponsive = false(size(data,2),1);
responses = struct;
tic
for iT = 1:size(data,2)
    
    % get info
    nCells = size(data(iT).allTrace,2)-1; % exclude recording cell
    nSweeps = size(data(iT).allTrace,3) * length(stimIdx);
    nPts = length(-preBuff:postBuff-1);    

    % initialize 
    exCellIdx = false(nCells,1);
    inCellIdx = false(nCells,1);
    nrCellIdx = false(nCells,1);
    % patchResp = zeros(nPts,nSweeps);
    cellResps = zeros(nPts,nSweeps,nCells);
    cellDist = zeros(nCells,1);

    % check patched cell
    cellTrace = squeeze(data(iT).allTrace(:,1,:));
    [patchResp,patchExcited,~,~] = test_cell_response(cellTrace,stimIdx,preBuff,postBuff);
    if patchExcited
        patchedResponsive(iT) = true;
    end

    % go through the rest of the cells 
    for iCell = 1:nCells
        
        % check response 
        cellTrace = squeeze(data(iT).allTrace(:,iCell+1,:)); % +1 to ignore the patched cell 
        [cellResps(:,:,iCell),exCellIdx(iCell),inCellIdx(iCell),nrCellIdx(iCell)] = ...
            test_cell_response(cellTrace,stimIdx,preBuff,postBuff);

        % scale distance 
        cellDist(iCell) = data(iT).dist(iCell+1) * pix2micron;  % +1 to ignore the patched cell 
    end 

    % output
    responses(iT).patchResp = patchResp;
    responses(iT).exCellIdx = exCellIdx;
    responses(iT).inCellIdx = inCellIdx;
    responses(iT).nrCellIdx = nrCellIdx;
    responses(iT).cellTraces = cellResps;
    responses(iT).cellDist = cellDist;
    responses(iT).perE = (sum(exCellIdx)) / nCells; 
    responses(iT).perI = (sum(inCellIdx)) / nCells; 
    responses(iT).perN = (sum(nrCellIdx)) / nCells; 
    responses(iT).distE = cellDist(exCellIdx);
    responses(iT).distI = cellDist(inCellIdx);
    responses(iT).distN = cellDist(inCellIdx);
end

% only include the recordings where the patched cell responded 
responses = responses(patchedResponsive);
toc
end 

%% 

function [snippet,excited,inhibited,noResponse] = test_cell_response(cellTrace,stimIdx,preBuff,postBuff)

% initialize 
excited = false; inhibited = false; noResponse = false; 
snippet = []; pre = []; post = [];

% go through each sweep and stim 
for iSweep = 1:size(cellTrace,2)
    for iStim = 1:length(stimIdx)
        tmp = cellTrace(stimIdx(iStim)-preBuff:stimIdx(iStim)+postBuff-1,iSweep);
        tmp = smoothdata(tmp,'gaussian',10);
        tmp = (tmp - mean(tmp(1:preBuff))) ./ std(tmp(1:preBuff));
        %         tmp = (tmp - mean(tmp))  ./ std(tmp);
        snippet = cat(2,snippet,tmp);
        pre = cat(1,pre,median(tmp(1:preBuff)));
        post = cat(1,post,median(tmp(preBuff+1:preBuff*2)));
    end
end

% compute significance
[p,h] =  signrank(pre,post);

% assign excitation vs inhibition
cm = mean(snippet,2); % mean trace
if h
    if median(cm(1:preBuff)) < median(cm(preBuff+1:preBuff*2))
        excited = true;
    else
        inhibited = true;
    end
else
    noResponse = true;
end
end