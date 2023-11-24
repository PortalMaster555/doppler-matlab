function fourier = pos_fft(inputRegion)
% POSFFT Computes fft, then finds true amplitudes for positive frequencies.
    N = size(inputRegion, 1);
    y = fft(inputRegion);
    absY = abs(y/N); %scaled two-sided amplitudes in f-space
    if mod(N,2) == 0
        fourier = absY(1:N/2+1); % take the second half
    else
        fourier = absY(1:(N+1)/2); % take the second half
    end
    fourier(2:end-1) = 2*fourier(2:end-1); % double all but 0 and nyq.
end