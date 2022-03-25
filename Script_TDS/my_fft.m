% v2.0 - now one can use the freq_min = 0, freq_max <=0 to keep the
% machine-decided size of FFT (but still only on the positive freq halfspace)
function [freq, ampl_f, phase_f] = my_fft(t, ampl_t, freq_min, freq_max, phase_unwrap_lim)


TimeSignalSizePOT=nextpow2(length(ampl_t));
TimeSignalSize=2^TimeSignalSizePOT;
% TimeSignalSize=2^8;

 
 FourierTransform = fft(ampl_t, TimeSignalSize);
%FourierTransform = fft(ampl_t);


[zzz, FFTSize]=size(FourierTransform);


dt=abs(t(1)-t(2));
freq = (1/dt/FFTSize).*(0: (FFTSize/2 - 1));

[zzz NewFFTSize] = size(freq); % only DC and positive freqs space is to be returned

ampl_f = abs(FourierTransform);
% phase_f = -unwrap(angle(FourierTransform), phase_unwrap_lim);
phase_f = -unwrap(angle(FourierTransform), phase_unwrap_lim);

ampl_f=ampl_f(1:NewFFTSize)./(TimeSignalSize); % divide by TimeSignalSize to get a correctly normalized spectral density
phase_f=phase_f(1:NewFFTSize);

% Forming the output within given freq boundaries
% MinFreqInd = find(freq >= freq_min, 1);
% MaxFreqInd = find(freq >= freq_max, 1);


if freq_min==0
    MinFreqInd = 1
else
MinFreqInd = find(freq >= freq_min, 1);
end

if freq_max <= 0
MaxFreqInd = NewFFTSize;
else
MaxFreqInd = find(freq >= freq_max, 1);
end

ampl_f=ampl_f(MinFreqInd:MaxFreqInd);
phase_f=phase_f(MinFreqInd:MaxFreqInd);
freq=freq(MinFreqInd:MaxFreqInd);


% deltaf=1/deltat;
% dt=1/HalfFreqSpaceSize/deltat;
% freq = deltaf*(0:FFTSize)/HalFreqSpaceSize;

% time_span=abs(max(t)-min(t));

% dfreq=1/time_span; % sampling frequency


% freq=dfreq.*(0: FFTsize/2)./FFTsize;


