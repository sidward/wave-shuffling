function [psf, py, pz] = wave_coeffs_to_psf(c, l, wx, y, z, scl)
  % [psf] = wave_coeffs_to_psf(c, l, wx, y, z, nodc)
  %
  % Convert wave coefficients to wave psf.
  %
  % Inputs:
  %   c    Sparse FFT coefficients that represent the sinusoidal gradient waveform. 
  %   l    Locations of the sparse FFT coefficients.
  %   wx   Size of the readout dimension after wave spreading.
  %   y    Spatial positions of the first phase encode dimension.
  %   z    Spatial positions of the second phase encode dimension.
  %   scl  Scale DC term.
  %
  % Outputs:
  %   psf  Wave point spread function.

  if (nargin < 6)
    scl = 1;
  end

  % Data dimensions
  sy = numel(y);
  sz = numel(z);

  % Initializing
  pY = zeros(wx, 1);
  pZ = zeros(wx, 1);

  % Extracting Y and Z positions.
  cY = c(1:end/2);
  cZ = c(1 + end/2:end);
  lY = l(1:end/2);
  lZ = l(1 + end/2:end);

  % Extracting real and imaginary parts.
  pY(lY) = cY(1:end/2) + 1i * cY(1 + end/2:end);
  pZ(lZ) = cZ(1:end/2) + 1i * cZ(1 + end/2:end);

  % Scaling DC.
  pY(lY == 1) = scl * pY(lY == 1);
  pZ(lZ == 1) = scl * pZ(lZ == 1);

  % Conjugate symmetry
  pY(mod(wx - lY + 1, wx) + 1) = conj(pY(lY));
  pZ(mod(wx - lZ + 1, wx) + 1) = conj(pZ(lZ));

  % Getting phase per cm.
  py = real(ifft(pY(:)));
  pz = real(ifft(pZ(:)));

  % Getting PSF arrays
  psfY = wave_psf(wx, y, py);
  psfZ = wave_psf(wx, z, pz); psfZ = reshape(psfZ, size(psfZ, 1), 1, size(psfZ, 2));
  psf  = bsxfun(@times, psfY, psfZ);
end
