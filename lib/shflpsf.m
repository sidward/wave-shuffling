function [w, c, shift, locs, res] = shflpsf(x, w, y_axis, z_axis, c_init, locs, num_iter, sigmoid_lambda, tol)
  % [w, c, shift, locs, res] = shflpsf(x, w, y_axis, z_axis, c_init, locs, num_iter, sigmoid_lambda, tol)
  %
  % Calibrate for hardware related wave-encoding errors as a deconvolution step.
  %
  % Inputs:
  %   x              Wave reconstruction using theoretical wave-PSF.
  %   w              Theoretical wave-PSF.
  %   y_axis         Spatial axis along the phase-encode direction.
  %   z_axis         Spatial axis along the partition-encode direction.
  %   c_init         Fourier coefficients of theoretical wave sinusoid.
  %   locs           Fourier locations of theoretical wave sinusoid.
  %   num_iter       Alternate between coefficient optimization and spatial offset optimization.
  %   sigmoid_lamda  Scaling term to determine edge detector sensitivity.
  %   tol            Optimization tolerance.
  %
  % Outputs:
  %   w              Corrected wave-PSF.
  %   c              Corrected Fourier coefficients.
  %   shift          Estimated spatial offset.
  %   locs           Locations of Fourier coefficients.
  %   res            Deconvolution result.

  shift = zeros(2, num_iter);
  c     = 1 * c_init;
  [w, ~, c, ~, ~, r]                          = globalpsf(x, w, y_axis, z_axis, c, locs, sigmoid_lambda, tol);
  for iter = 1:num_iter
    [w, ~, shift(1, iter), shift(2, iter), r] = offsetpsf(r, w, y_axis, z_axis, c, locs, sigmoid_lambda, tol);
    [w, ~, c, ~, ~, r]                        = globalpsf(r, w, y_axis, z_axis, c, locs, sigmoid_lambda, tol);
  end
  res = r;
end
