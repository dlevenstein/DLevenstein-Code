function [ pupildilation ] = GetPupilDilation( baseName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Note: you may need gstreamer: https://gstreamer.freedesktop.org/download/

baseName = 'Layers_LFP_Test02_170323_151411';
recordingsfolder = 'C:\Users\rudylabadmin\Desktop\layers\Recordings';
recordingsfolder = '/mnt/proraidDL/Database/WMProbeData';

addpath(fullfile(recordingsfolder,baseName));

vidName = fullfile(recordingsfolder,baseName,[baseName,'.avi']);
abfname = fullfile(recordingsfolder,baseName,[baseName,'.abf']);
analogName = fullfile(recordingsfolder,baseName,['analogin.dat']);

savefile = fullfile(recordingsfolder,baseName,[baseName,'.pupildiameter.behavior.mat']);
savevid = fullfile(recordingsfolder,baseName,[baseName,'.pupilvid.avi']);

SAVEVID = true;
savevidfr = 10;
if SAVEVID
    pupdiamVid = VideoWriter(savevid);
    pupdiamVid.FrameRate = 1./(0.015.*savevidfr);
    open(pupdiamVid);
end
%%
%Load in the video (has to be in not .avi for now...)
pupilvidobj = VideoReader(vidName);
NumberOfFrames = pupilvidobj.NumberOfFrames;

%Use the first frame to get image dimensions
vidframe=read(pupilvidobj,10);
imagesize=size(vidframe);

%Initiate the pupil area vector, center coordinates
puparea = nan(NumberOfFrames,1);
pupcoords = nan(NumberOfFrames,2);

%Loop through the frames and get pupil location/area
 figure
 colormap('gray')
for ff = 1:NumberOfFrames 
%ff = 1;
if mod(ff,200)==0
display([num2str(ff),' of ' num2str(NumberOfFrames),' Frames'])
end

     vidframe=read(pupilvidobj,ff);
     
     %Convert to greyscale
    if size(vidframe,3)==3
        vidframe = rgb2gray(vidframe);
    end
    vidframe_orig = vidframe;


    
    %Define the mask: trace the eye
    if ff==1
        
        subplot(2,2,1)
        imagesc(vidframe_orig);
        title('Trace the Eye')
        h = imfreehand;
        mask = ~createMask(h);
        %keyboard
%         mask = false(size(vidframe));
%         ymax = 70; ymin = 40; xmax = 60; xmin = 35;
%         mask(xmax:end,:) = true;
%         mask(1:xmin,:) = true;
%         mask(:,1:ymin) = true;
%         mask(:,ymax:end) = true;
%         
%         %
%         maskval =  max(vidframe(:));
        maskval = 255;
    end
    vidframe(mask) = maskval;
    

    %Define the pupil size and intensity threshold - trace the pupil
    if ff==1
        
            %Show the Masked frame
            subplot(2,2,2)
         imagesc(vidframe)
        title('Trace the Pupil')
        
        h = imfreehand;
        pupilmask = createMask(h);
        pupilpixels = vidframe(pupilmask); 
        irispixels  = vidframe(~mask & ~ pupilmask); 
        
       % keyboard
        pupilsizethresh = 50; %Pupil must be larger than this many pixels.
        
        % Histograms of iris and pupil
%         figure
%             hist(double(irispixels)./255)
%             hold on
%             hist(double(pupilpixels)./255)
%             plot(mean(double(pupilpixels)./255).*[1 1],get(gca,'ylim'),'k')
%             plot(mean(double(pupilpixels)./255)+1.*std(double(pupilpixels)./255).*[1 1],get(gca,'ylim'),'k')
        %%
        %Pupil must be darker than this intensity: 2.5std above mean pupil
        intensitythresh = mean(double(pupilpixels)./255)+2.*std(double(pupilpixels)./255); 
    end
        
    %Get the pupil - all pixels below intensity threshold (black)
    pupil=~im2bw(vidframe,intensitythresh);

    if SAVEVID && mod(ff,savevidfr)==0; 
        subplot(4,4,9);imagesc(pupil); end
    %Show the thresholded image

    
    %Remove small objects and holes (i.e. dark/bright spots)
    pupil=bwmorph(pupil,'close');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,10);imagesc(pupil); end
    
    pupil=bwmorph(pupil,'open');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,11);imagesc(pupil); end
    
    pupil=bwareaopen(pupil,pupilsizethresh);
    pupil=imfill(pupil,'holes');
    if SAVEVID && mod(ff,savevidfr)==0; 
         subplot(4,4,12);imagesc(pupil); end

    % Tagged objects in BW image
    L=bwlabel(pupil);
    % Get areas and tracking rectangle
    out_a=regionprops(L);
    % Count the number of objects
    N=size(out_a,1);
    
    if N < 1 || isempty(out_a) % Returns if no object in the image
        continue
    end

    % Select larger area
    areas=[out_a.Area];
    [area_max pam]=max(areas);
    puparea(ff) = area_max;
    
    centro=round(out_a(pam).Centroid);
    X=centro(1);
    Y=centro(2);
  
    if SAVEVID && mod(ff,savevidfr)==0;  
    subplot(2,2,1)
    imagesc(vidframe_orig);
    hold on
    rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],...
        'Curvature', [1,1],'LineWidth',1)
    plot(X,Y,'g+')
    hold off
    end
    
    pupcoords(ff,1) = X; pupcoords(ff,2) = Y;
    
   if SAVEVID && mod(ff,savevidfr)==0;  
   subplot(4,1,4)
   plot(1:ff,puparea(1:ff),'k')
   end
     
   if SAVEVID && mod(ff,savevidfr)==0;
       imgFrame = getframe(gcf);
       writeVideo(pupdiamVid,imgFrame.cdata);
   end
end

close(pupdiamVid)

%0-1 Normalize
puparea_pxl = puparea;
puparea = (puparea-min(puparea))./(max(puparea)-min(puparea));

%% Load the .abf for the timestamps

timepulses = readmulti(analogName,1);

sf_pulse = 1./20000; %Sampling Frequency of the .abf file
t_pulse = [1:length(timepulses)]'.*sf_pulse;

pulsethreshold =1e4;  %Adjust this later to set based on input.
pulseonsets = find(diff(timepulses<pulsethreshold)==1);
pulset = t_pulse(pulseonsets);

expectedpulserate = 0.015;
%%
shortpulses=diff(pulset)<(0.5.*expectedpulserate);
pulset(shortpulses) = [];

% longpulses=diff(pulset)>(2.*expectedpulserate);
% pulset(longpulses) = [];

%hist(diff(pulset))

%%
% timepulses = abfload(abfname);
% timepulses = timepulses(:,1);
% 
% sf_abf = 1./20000; %Sampling Frequency of the .abf file
% t_abf = [1:length(timepulses)]'.*sf_abf;
% 
% pulsethreshold =1;  %Adjust this later to set based on input.
% pulseonsets = find(diff(timepulses<pulsethreshold)==1);
% pulset = t_abf(pulseonsets);
% 
% 
% pulset(1) = []; %remove the first trigger... make this more rigorous later 

%Check that the number of identified pulses = the number of frames
if length(pulset)~=length(puparea); 
    display('WARNING: NUMBER OF FRAMES DON"T MATCH!!!    v sad.');%keyboard; 
end

%%
figure
plot(t_pulse,timepulses,'k')
hold on
plot(pulset,zeros(size(pulset)),'r+')

%%
figure
plot(puparea)
%%





%% Behavior struct

pupildilation.t_abf = pulset;
%pupildilation.t_interp;
pupildilation.puparea = puparea;
pupildilation.puparea_pxl = puparea_pxl;
pupildilation.pupcoords = pupcoords;
pupildilation.mask = mask;


save(savefile,'pupildilation')


end

