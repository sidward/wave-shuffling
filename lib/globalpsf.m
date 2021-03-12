function [w_new, w_err, c_new, c_err, l, res] = globalpsf(x, w, y_axis, z_axis, c_init, fft_locations, sigmoid_lambda, tol)
  % Please see shflpsf for documentation.

  % Helper functions.
  fl  = @(x) x(:);
  nrm = @(x) x/norm(x(:));

  % Dimensions.
  wx = size(w, 1);
  sx = size(x, 1);
  sy = size(x, 2);
  sz = size(x, 3);
  tk = size(x, 4);

  % Unpack locations.
  l = [fl(fft_locations(1:end/2)); fl(fft_locations(1 + end/2:end))];

  % Initialize coefficients.
  c = [fl(zeros(2 * numel(l), 1))];

  % Get isocenter.
  [~, center_y] = sort(abs(y_axis), 'ascend');
  [~, center_z] = sort(abs(z_axis), 'ascend');
  center_y = center_y(1);
  center_z = center_z(1);

  % Extract slices and reshape.
  y = zeros(wx, sy, sz, tk);
  y((wx - sx)/2 + [1:sx], :, :, :) = x;
  x = cat(4, squeeze(y(:, :, center_z, :)), squeeze(y(:, center_y, :, :)));

  % Normalize.
  for p = 1:size(x, 3)
    for q = 1:size(x, 4)
      x(:, :, p, q) = nrm(x(:, :, p, q));
    end
  end

  % Define objective function.
  function j = objective(c);
    % Unpack Y and Z coeffs.
    cY = c * 1; cY(1 + end/2:end  ) = 0;
    cZ = c * 1; cZ(        1:end/2) = 0;

    % Estimate PSFs.
    wY = squeeze(wave_coeffs_to_psf(cY, l, wx, y_axis, 0));
    wZ = squeeze(wave_coeffs_to_psf(cZ, l, wx, 0, z_axis));

    % Prepare error PSFs.
    we = cat(4, wY, wZ);

    % Edge image.
    prj = F_inv(bsxfun(@times, we, F_fwd(x, 1)), 1);
    dlt = prj - circshift(prj, 1, 1);

    % Amplify edges.
    sigmoid = @(x) (1 ./ (1 + exp(-x))) - 0.5;
    if (sigmoid_lambda > 0)
      dlt = sigmoid(sigmoid_lambda * abs(dlt));
    end

    % Cost function.
    j = norm(fl(dlt), 1);
  end

  % Perform optimization.
  opt = [];
  opt.Display     = 'iter';
  opt.MaxIter     = 100000;
  opt.MaxFunEvals = 100000;
  opt.TolX        = tol;
  opt.TolFun      = tol;
  fprintf('---------------------------------------------------------------------------------------------\n');
  fprintf('> GLOBAL PSF.\n');
  fprintf('---------------------------------------------------------------------------------------------\n');
  c_err = fminsearch(@(c) objective(c), c, opt);
  c_new = c_init(:) - c_err(:);

  % Estimate error PSFs.
  w_err = wave_coeffs_to_psf(c_err, l, wx, y_axis, z_axis);
  w_new = wave_coeffs_to_psf(c_new, l, wx, y_axis, z_axis);

  % Prepare result.
  y   = F_inv(bsxfun(@times, w_err, F_fwd(y, 1)), 1);
  y   = y((wx - sx)/2 + [1:sx], :, :, :);
  res = y/norm(y(:));
end
