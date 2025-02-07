function figHandle = init_figure(widthCm, heightCm, varargin)
% init_figure  Create a figure with specified width/height in cm, and set font defaults.
%
%   figHandle = init_figure(widthCm, heightCm)
%   figHandle = init_figure(widthCm, heightCm, 'Name', Value, ...)
%
% REQUIRED INPUTS:
%   widthCm  - figure width (in centimeters)
%   heightCm - figure height (in centimeters)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'FontName' (char) : Font for axes labels and tick labels (default: 'Arial')
%   'FontSize' (numeric): Font size for axes labels and tick labels (default: 14)
%
% PROPERTIES SET AUTOMATICALLY:
%   Figure Units            : 'centimeters'
%   PaperPositionMode       : 'auto'
%   Figure Background Color : 'w'
%
% EXAMPLE:
%   % Create a 12 cm x 8 cm figure, using default Arial 14
%   f = initFigure(12, 8);
%
%   % Create a 10 cm x 6 cm figure with Times New Roman, font size 16
%   f = initFigure(10, 6, 'FontName', 'Times New Roman', 'FontSize', 16);
%
%   % Plot something
%   ax = axes('Parent', f);
%   plot(ax, 1:10, rand(1,10), 'o-');
%   xlabel(ax, 'X Label');
%   ylabel(ax, 'Y Label');
%   title(ax, 'Example Plot');

    %------------------------------
    % 1) Parse inputs
    %------------------------------
    parser = inputParser();
    addRequired(parser, 'widthCm',  @(x) isnumeric(x) && isscalar(x) && x>0);
    addRequired(parser, 'heightCm', @(x) isnumeric(x) && isscalar(x) && x>0);
    
    addParameter(parser, 'FontName', 'Arial', @ischar);
    addParameter(parser, 'FontSize', 14,      @(x) isnumeric(x) && isscalar(x));
    
    parse(parser, widthCm, heightCm, varargin{:});
    
    widthCm   = parser.Results.widthCm;
    heightCm  = parser.Results.heightCm;
    fontName  = parser.Results.FontName;
    fontSize  = parser.Results.FontSize;

    %------------------------------
    % 2) Create the figure
    %------------------------------
    figHandle = figure('Units',               'centimeters', ...
                       'Position',            [2, 2, widthCm, heightCm], ...
                       'PaperPositionMode',   'auto', ...
                       'Color',               'w');
                   
    %------------------------------
    % 3) Set default font properties
    %    for any axes/text created in this figure
    %------------------------------
    set(figHandle, 'DefaultAxesFontName',  fontName, ...
                   'DefaultAxesFontSize',  fontSize, ...
                   'DefaultTextFontName',  fontName, ...
                   'DefaultTextFontSize',  fontSize);

   hold on
end
