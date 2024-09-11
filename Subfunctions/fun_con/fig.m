function [h_fig, wid_mon, hei_mon] = fig
% PURPOSE: To generate a fullscreen figure
%---------------------------------------------------
% USAGE: fig
% --------------------------------------------------
% SEE ALSO: f(results)
%---------------------------------------------------
% REFERENCES: 
%---------------------------------------------------
% REMARKS: 
% h.PaperUnits = 'centimeters';
% (32.4/1200)*1080
% h.PaperSize = [51.84 29.16];
% 51.84 cm x 32.4 cm
%---------------------------------------------------

% Written by: Kadir Mutlu

wid_mon = 1920;
hei_mon = 1080;
h_fig = figure('Color', 'white', 'Position', [0 0 wid_mon hei_mon]);
h_fig.PaperUnits = 'centimeters';
h_fig.PaperSize = [51.84 29.16];% real size of the fig.
h_fig.Renderer = 'painters';
end
