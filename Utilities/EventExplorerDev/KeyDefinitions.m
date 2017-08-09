function KeyDefinitions(f, e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Get the gui data
obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj); 

switch e.Key
    case 'uparrow'
        FO.scaleLFP = FO.scaleLFP+0.1;
        EventVewPlot
    case 'downarrow'
        FO.scaleLFP = FO.scaleLFP-0.1;
        EventVewPlot
end

guidata(FO.fig, FO); %write any new data to the object
end

