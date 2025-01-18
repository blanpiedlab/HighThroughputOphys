
function yout = gaussfilter(yin, sampfreq, cornerfreq)
% YOUT = GAUSSFILTER(YIN, SAMPFREQ, CORNERFREQ)
% returns the Gaussian filtered version of YIN. SAMPFREQ is the sampling frequency
% CORNERFREQ is the corner frequency of the Gaussian
% HPCR, 1999
% Example usage: filter a sweep acquired at 20000 Herz at a corner frequency of 1 kHz:
% filtered = gaussfilter(sweep, 20000, 1000)
fc = cornerfreq/sampfreq;  % corner frequency in units of the sampling frequency
sigma = 0.132505/fc;
nc = round(4*sigma); % number of coefficients on each side, not counting central one
coeffs = -nc:nc;
coeffs = exp((-coeffs.^2)/(2*sigma^2))/(sqrt(2*pi)*sigma);
yout = conv(yin,coeffs);
yout = yout(nc+1:end-nc);  % chop off ends so that same length as yin
% changed as a result of Ingo's problem  9/1/01
% yout = convn(yin,coeffs,'same');   % seems to be a bug in convn!
end