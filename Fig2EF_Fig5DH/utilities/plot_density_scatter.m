% function plot_density_scatter(data, xOffset, color)
%     [f, ~] = ksdensity(sort(data));
%     x = linspace(0, length(f), length(f));
%     f = rescale(f, 0, 0.25);
%     xi = linspace(0, length(f), length(data));
%     f = interp1(x, f, xi);
%     x = rand(1, length(data));
%     x = rescale(x, 0, 1);
%     x = x .* f;
%     x(2:2:end) = -x(2:2:end);
%     x = x + xOffset;
%     scatter(x, sort(data), 20, color, 'filled');
% end

function plot_density_scatter(data, xOffset, color, arrangement)
    % Compute density estimation
    [f, xi] = ksdensity(sort(data));
    f = rescale(f, 0, 0.25);
    
    % Interpolate density for the actual data points
    dataSorted = sort(data);
    f_interp = interp1(xi, f, dataSorted, 'linear', 'extrap');
    
    % Structured placement
    n = length(data);
    x = zeros(1, n);
    
    switch arrangement
        case 'up'
            x = (1:n) / n .* f_interp;  % Sequentially spread upwards
        case 'down'
            x = -(1:n) / n .* f_interp; % Sequentially spread downwards
        case 'fan'
            x(1:2:end) = (1:ceil(n/2)) / ceil(n/2) .* f_interp(1:2:end);
            x(2:2:end) = -(1:floor(n/2)) / floor(n/2) .* f_interp(2:2:end);
        case 'square'
            numRows = ceil(sqrt(n));
            rowIndex = mod(0:n-1, numRows) + 1;
            x = (rowIndex - median(rowIndex)) / numRows .* max(f_interp);
        case 'hex'
            numRows = ceil(sqrt(n));
            rowIndex = mod(0:n-1, numRows) + 1;
            shift = mod(0:n-1, 2) * 0.5; % Hex pattern shift
            x = (rowIndex - median(rowIndex)) / numRows .* max(f_interp) + shift;
        case 'rand'  % Preserve the original random layout
            x = rand(1, n) .* f_interp;
            x(2:2:end) = -x(2:2:end);
        otherwise
            error('Invalid arrangement option');
    end
    
    x = x + xOffset;
    scatter(x, dataSorted, 20, color, 'filled');
end
