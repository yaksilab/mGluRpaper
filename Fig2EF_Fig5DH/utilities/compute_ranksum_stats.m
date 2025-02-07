function [stats] = compute_ranksum_stats(dataSet1, dataSet2, statName)
% compute_ranksum_stats - Computes rank-sum statistics and additional information.
%
% Inputs:
%   dataSet1   - Group 1 data
%   dataSet2   - Group 2 data
%   statName   - A name or identifier for the stats structure
%
% Outputs:
%   stats      - A structure containing all relevant test statistics
%
% Example Usage:
%   stats = compute_ranksum_stats(group1, group2, 'Comparison');

    % Define a helper function for standard error of the mean (SEM)
    SEM = @(x) (nanstd(x) / sqrt(length(x)));

    % Calculate means and SEMs for both groups
    avgGroup1 = round(nanmean(dataSet1), 3, "significant");
    semGroup1 = round(SEM(dataSet1), 2, "significant");

    avgGroup2 = round(nanmean(dataSet2), 3, "significant");
    semGroup2 = round(SEM(dataSet2), 2, "significant");

    % Perform rank-sum test (Mann-Whitney U test)
    [p, h, statsInfo] = ranksum(dataSet1, dataSet2);

    % Extract statistics if available
    zval = NaN;
    rankSum = NaN;
    
    if ~isempty(statsInfo)
        if isfield(statsInfo, 'zval')
            zval = statsInfo.zval;
        end
        if isfield(statsInfo, 'ranksum')
            rankSum = statsInfo.ranksum;
        end
    end

    % Create the output structure
    stats = struct( ...
        'name', statName, ...
        'p', p, ...                  % p-value
        'h', h, ...                  % Test decision (0 or 1)
        'zval', zval, ...            % Z-statistic
        'rankSum', rankSum, ...      % Rank sum
        'group1n', length(dataSet1), ...  % Sample size for group 1
        'group1', dataSet1, ...           % Raw data for group 1
        'group1Avg', avgGroup1, ...       % Mean for group 1
        'group1Sem', semGroup1, ...       % SEM for group 1
        'group2n', length(dataSet2), ...  % Sample size for group 2
        'group2', dataSet2, ...           % Raw data for group 2
        'group2Avg', avgGroup2, ...       % Mean for group 2
        'group2Sem', semGroup2 ...        % SEM for group 2
    );
end
