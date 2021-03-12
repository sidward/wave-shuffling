function y = F_fwd(x, ax)
  % y = F_fwd(x, ax) 
  %
  % Forward centered unitary Fourier transform along axes in ax.
  %
  % Inputs:
  %   x   Input array.
  %   ax  Axes to take the Fourier transform of.
  %
  % Output:
  %   y   Result.

  y = fftshift(fft(ifftshift(x, ax(1)),  size(x, ax(1)), ax(1)), ax(1)) ./ sqrt(size(x, ax(1)));
  for k = 2:numel(ax)
    y = fftshift(fft(ifftshift(y, ax(k)),  size(y, ax(k)), ax(k)), ax(k)) ./ sqrt(size(y, ax(k)));
  end

end
