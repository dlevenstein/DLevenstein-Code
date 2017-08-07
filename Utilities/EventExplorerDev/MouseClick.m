function MouseClick(e, src)
    obj = findobj('tag','EventExplorerMaster');  FO = guidata(obj); 
    
    %Get the info about the click with respect to figure, axes
    mousepoint = get(gcf, 'CurrentPoint');
    clicktype = get(gcf, 'SelectionType');
    clickcoords = get(gca, 'CurrentPoint');clickcoords = clickcoords(1,1:2);
    
    %Determine what the click does as a function of current action
    switch FO.currentuseraction
        
        case 'MarkEvents'
            switch clicktype
                case 'normal' %left mouse click
                    plot(clickcoords(1),clickcoords(2),'ro')
                    %Add the event to the list of marked events 1=leftclick
                    FO.markedevents = [FO.markedevents;...
                        [clickcoords 1]];
                case 'alt' %right mouse click
                    plot(clickcoords(1),clickcoords(2),'rx')
                    %Add the event to the list of marked events 3=rghtclick
                    FO.markedevents = [FO.markedevents;...
                        [clickcoords 3]];
            end
            
                      
    end
    
    guidata(FO.fig, FO); %write any new data to the object
end