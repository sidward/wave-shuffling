function psf = wave_psf(wx, locs, ppcm)
  % psf = wave_psf(wx, locs, ppcm)
  %
  % Generates wave PSF to use in a forward model.
  %
  % Inputs:
  %   wx    Readout dimension after accounting for wave-based aliasing.
  %   locs  Axis locations in cm.
  %   ppcm  See phase_per_cm_cosine..
  %
  % Output:
  %   psf   Wave PSF in hybrid space.
  psf = zeros(wx, numel(locs));
  for idx=1:numel(locs)
    psf(:, idx) = exp(-1i .* ppcm .* locs(idx));
  end
  psf = psf./abs(psf);
end
