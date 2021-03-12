function [w_new, w_err, shift_y, shift_z, res] = offsetpsf(x, w, y_axis, z_axis, c_init, l, sigmoid_lambda, tol)
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

  % Initialize coefficients.
  c = zeros(2, 1);

  % Get isocenter.
  [~, center_y] = sort(abs(y_axis), 'ascend');
  [~, center_z] = sort(abs(z_axis), 'ascend');
  center_y = center_y(1);
  center_z = center_z(1);

  % Extract PSF coefficients.
  cY = c_init .* 1; cY(1+end/2:end  ) = 0;
  cZ = c_init .* 1; cZ(      1:end/2) = 0;

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
    % Unpack.
    shift_y = c(1);
    shift_z = c(2);

    % Estimate PSFs.
    wY = squeeze(wave_coeffs_to_psf(cY, l, wx, shift_y * ones(size(y_axis)), 0));
    wZ = squeeze(wave_coeffs_to_psf(cZ, l, wx, 0, shift_z * ones(size(z_axis))));

    % Prepare error PSFs.
    we = cat(4, wY, wZ);

    % Edge image.
    prj = F_inv(bsxfun(@times, we, F_fwd(x, 1)), 1);
    dlt = prj - circshift(prj, 1, 1);

    % TODO: TEST
    dlt = dlt - circshift(dlt, 1, 2);

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
  fprintf('> OFFSET PSF.\n');
  fprintf('---------------------------------------------------------------------------------------------\n');
  c = fminsearch(@(c) objective(c), c, opt);
  shift_y = c(1);
  shift_z = c(2);
  fprintf('---------------------------------------------------------------------------------------------\n');
  fprintf('Shift Y: %f cm\n', shift_y);
  fprintf('Shift Z: %f cm\n', shift_z);

  % Estimate PSFs.
  wY = reshape(squeeze(wave_coeffs_to_psf(cY, l, wx, shift_y * ones(size(y_axis)), 0)), wx, sy, 1);
  wZ = reshape(squeeze(wave_coeffs_to_psf(cZ, l, wx, 0, shift_z * ones(size(z_axis)))), wx, 1, sz);

  w_err = wY .* wZ;
  w_new = w .* conj(w_err);

  % Prepare result.
  y   = F_inv(bsxfun(@times, w_err, F_fwd(y, 1)), 1);
  y   = y((wx - sx)/2 + [1:sx], :, :, :);
  res = y/norm(y(:));
  fprintf('---------------------------------------------------------------------------------------------\n');
end
