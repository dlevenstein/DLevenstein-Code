function [ ] = ColorbarWithAxis(cbounds,label)
%ColorbarWithAxis(cbounds,label) adds a colorbar to a imageplot and labels
%axis and bounds with > and <.
%
%TO DO
%   change inputs to 'parms' with inputParser
%   add 3std functionality
%
%DLevenstein 2016
%% inputParse for input options




%%
cb = colorbar;

if isstring(cbounds); switch cbounds
	case '3std'
end; end
            

caxis(cbounds);
title(cb, label);

cb.Ticks = [cbounds(1) mean(cbounds) cbounds(2)];
cb.TickLabels = {['> ',num2str(cbounds(1))],mean(cbounds),...
    ['< ',num2str(cbounds(2))]};

end

