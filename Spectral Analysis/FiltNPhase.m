function [ filt_data, power, phase ] = FiltNPhase(data,fbounds,sampfreq,varargin)
%[filt_data,power,phase] = FiltNPhase(data,fbounds,sampfreq)
%
%Optional Arguement: filter order (num_cyc)
%
%TO DO:
%   -Add option for high/low pass filter instead of bandpass
%   -Improve input parsing for optional inputs
%
%Last Updated: 2/22/16
%DLevenstein
%%


%Filter order to pass to fir1. higher order improves frequency resolution
%but also increases computing time substantially.
%default: 3 cycles of low bound freq
if nargin == 4
    num_cyc = varargin{1};
else
    num_cyc = 3;
end

low_bound = 1/fbounds(1);   %Low bound period
filt_order = num_cyc*low_bound*sampfreq;    
N = ceil(filt_order);       %order must be integer

%Make bandpass FIR filter
Wn = fbounds./(0.5*sampfreq);    %Normalize frequency bounds to 1/2 sample 
                                %rate (nyquist frequency) for fir1.
B = fir1(N,Wn);                 %Designs an N'th order lowpass FIR digital
                                %filter and returns the filter coefficients
                                %in length N+1 vector B.

%Filter data and get phase using Hilbert transform
filt_data = filtfilt(B,1,data);
power = abs(hilbert(filt_data));
phase = angle(hilbert(filt_data));

end

