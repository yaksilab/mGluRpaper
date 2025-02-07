function add_significance_line(figHandle, ax, xStart, xEnd, pValue, varargin)
%ADD_SIGNIFICANCE_LINE creates a line that marks the statistical significance
%
% Draws a line that marks the statistical significance on the axes and
% returns the handles to the lines and text. This function uses the
% p-value thresholds:
%    - p > 0.05       -> no marker (line not drawn)
%    - 0.05 >= p > 0.01  -> '*'
%    - 0.01 >= p > 0.001 -> '**'
%    - p <= 0.001     -> '***'
%
% INPUTS
%   figHandle - handle to the figure
%   ax        - handle to the axes
%   xStart    - the start position (x) of the significance line
%   xEnd      - the end position (x) of the significance line
%   pValue    - the p-value used to decide the marker
%
% NAME-VALUE PAIRS (optional):
%   'height'       : Y position at which the line is drawn
%                    (if not provided, will be set automatically
%                     near the top of the current axis limits)
%
%   'marker'       : manually override the text shown on the line
%                    (default based on pValue thresholds)
%
%   'edgeLength'   : length of the small vertical ticks at xStart, xEnd
%                    (0 => no ticks). Default is 10% of (xEnd - xStart)
%
%   'markerSpace'  : vertical space between the line and the marker
%                    (default = 0.02 * yRange if not specified)
%
%   'color'        : color of the line and text (default = 'black')
%
%   'lineWidth'    : line width of the significance line (default = 1)
%
%   'fontWeight'   : 'bold' or 'normal' (default = 'normal')
%
%   'fontSize'     : font size of the text (default = 10)
%
%   'fontName'     : font name of the text (default = 'Arial')
%
% OUTPUT
%
%
% EXAMPLE:
%   figure; ax = axes; hold on;
%   x = 1:5; y = rand(1,5);
%   plot(ax, x, y, 'o-');
%   add_significance_line(gcf, ax, 1, 2, 0.003, 'color','r');
%

    %------------------------------
    % 1) Switch to figure / axes
    %------------------------------
    if ~isempty(figHandle) && ishandle(figHandle)
        figure(figHandle);  % bring that figure forward
    end
    axes(ax);
    hold(ax, 'on');
    
    %------------------------------
    % 2) Determine the default marker from p-value
    %------------------------------
    if pValue <= 0.001
        defaultMarker = '***';
    elseif pValue <= 0.01
        defaultMarker = '**';
    elseif pValue <= 0.05
        defaultMarker = '*';
    else
        % pValue > 0.05 -> no marker -> no line
        defaultMarker = '';
    end

    %------------------------------
    % 3) Parse optional inputs
    %------------------------------
    p = inputParser;
    addParameter(p, 'marker',      defaultMarker);
    addParameter(p, 'height',      []);       % automatically determined if empty
    addParameter(p, 'edgeLength',  []);       % default is 10% of (xEnd - xStart)
    addParameter(p, 'markerSpace', []);       % gap between line and marker
    addParameter(p, 'color',       'black');
    addParameter(p, 'lineWidth',   2);
    addParameter(p, 'fontWeight',  'normal');
    addParameter(p, 'fontSize',    20);
    addParameter(p, 'fontName',    'Arial');
    parse(p, varargin{:});
    
    marker      = p.Results.marker;
    lineHeight  = p.Results.height;
    edgeLength  = p.Results.edgeLength;
    markerSpace = p.Results.markerSpace;
    lineColor   = p.Results.color;
    lineWidth   = p.Results.lineWidth;
    fontWeight  = p.Results.fontWeight;
    fontSize    = p.Results.fontSize;
    fontName    = p.Results.fontName;

    % If no marker (i.e., pValue>0.05), return with empty handle
    if isempty(marker)
        return
    end
    
    % Default edgeLength to 10% of (xEnd - xStart) if not provided
    if isempty(edgeLength)
       edgeLength = 0.1 * (xEnd - xStart);
    end
    
    %------------------------------
    % 4) Auto-compute line height
    %    (if user didn't specify)
    %------------------------------
    yLim = get(ax, 'YLim');
    yRange = yLim(2) - yLim(1);
    if isempty(lineHeight)
        % Place the significance line near the top 10% of the axis range
        % but not off the plot
        offsetTop  = 0.1 * yRange;
        lineHeight = yLim(2) - offsetTop;
    end
    
    % If markerSpace was not specified, pick a small fraction of the y-range
    if isempty(markerSpace)
        markerSpace = 0.01 * yRange;  
    end
    
    %------------------------------
    % 5) Draw the main horizontal line
    %------------------------------
    hLine = plot(ax, [xStart, xEnd], [lineHeight, lineHeight], '-', ...
        'Color', lineColor, ...
        'LineWidth', lineWidth);

    %------------------------------
    % 6) Draw optional "edge" lines
    %------------------------------
    if edgeLength ~= 0
        hLeftEdge  = plot(ax, [xStart, xStart], [lineHeight-edgeLength, lineHeight], '-', ...
            'Color', lineColor, 'LineWidth', lineWidth);
        hRightEdge = plot(ax, [xEnd,   xEnd  ], [lineHeight-edgeLength, lineHeight], '-', ...
            'Color', lineColor, 'LineWidth', lineWidth);
    end
    
    %------------------------------
    % 7) Draw the marker (asterisks)
    %------------------------------
    hText = text(ax, 0.5*(xStart + xEnd), lineHeight + markerSpace, marker, ...
        'FontName',          fontName, ...
        'FontSize',          fontSize, ...
        'FontWeight',        fontWeight, ...
        'Color',             lineColor, ...
        'HorizontalAlignment','center');

end

