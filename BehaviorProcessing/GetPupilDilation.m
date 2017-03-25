function [ pupildilation ] = GetPupilDilation( baseName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Note: you may need gstreamer: https://gstreamer.freedesktop.org/download/

%baseName = '20170303_1';
recordingsfolder = 'C:\Users\rudylabadmin\Desktop\layers\Recordings';

addpath(fullfile(recordingsfolder,baseName));

vidName = fullfile(recordingsfolder,baseName,[baseName,'.avi']);
abfname = fullfile(recordingsfolder,baseName,[baseName,'.abf']);

savefile = fullfile(recordingsfolder,baseName,[baseName,'.pupildiameter.behavior.mat']);
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

    %Show the thresholded image
   % subplot(4,4,9)
   % imagesc(pupil)
    
    %Remove small objects and holes (i.e. dark/bright spots)
    pupil=bwmorph(pupil,'close');
    %subplot(4,4,10);imagesc(pupil)
    
    pupil=bwmorph(pupil,'open');
   % subplot(4,4,11);imagesc(pupil)
    
    pupil=bwareaopen(pupil,pupilsizethresh);
    pupil=imfill(pupil,'holes');
   % subplot(4,4,12);imagesc(pupil);

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
    
%     subplot(2,2,1)
%     imagesc(vidframe_orig);
%     hold on
%     rectangle('Position',out_a(pam).BoundingBox,'EdgeColor',[1 0 0],...
%         'Curvature', [1,1],'LineWidth',1)
%     plot(X,Y,'g+')
    
    pupcoords(ff,1) = X; pupcoords(ff,2) = Y;
    
    
  %  subplot(4,1,4)
  %  plot(1:ff,puparea(1:ff),'k')
     
end

%0-1 Normalize
puparea = (puparea-min(puparea))./(max(puparea)-min(puparea));

%% Load the .abf for the timestamps

timepulses = abfload(abfname);
timepulses = timepulses(:,1);

sf_abf = 1./20000; %Sampling Frequency of the .abf file
t_abf = [1:length(timepulses)]'.*sf_abf;

pulsethreshold =1;  %Adjust later to set based on input.
pulseonsets = find(diff(timepulses<pulsethreshold)==1);
pulset = t_abf(pulseonsets);

%Check that the number of identified pulses = the number of frames
if length(pulset)~=length(puparea); keyboard; end

%%
figure
plot(t_abf,timepulses,'k')
hold on
plot(pulset,zeros(size(pulset)),'r+')


%% Behavior struct

pupildilation.t = pulset;
pupildilation.puparea = puparea;
pupildilation.pupcoords = pupcoords;
pupildilation.mask = mask;


save(savefile,'pupildilation')


end

