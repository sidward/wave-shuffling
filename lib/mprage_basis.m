function [u, sig, sv, t1_vals] = mprage_basis(tk, t1_vals, ttfa, tr, alpha, esp, tf, reps, bias)
  % [u, sig, sv, t1_vals] = genBasis(tk, t1_vals, ttfa, tr, alpha, esp, tf, bias)
  %
  % Generates the shuffling-temporal basis from MPRAGE signal evolution simulations.
  %
  % Inputs:
  %   tk       Shuffling temporal dimension.
  %   t1_vals  T1 image map of dimensions (dimX, dimY, dimZ)
  %   ttfa     Time from 180* to first alpha in milliseconds.
  %   tr       Repitition time in milliseconds.
  %   alpha    Small tip flip angle in degrees.
  %   esp      Echo spacing, or time between adjacent alphas 
  %            in milliseconds.
  %   tf       Turbo factor, or the number of echos.
  %   reps     Number of times to repeat simulation to account for signal with large T1.
  %   bias     If true, include a basis vector to fit for bias.
  %
  % Outputs:
  %   u        Shuffling temporal basis.
  %   sig      Signal over time used to estimate the basis.
  %   sv       Singular values of the basis vectors.
  %   t1_vals  T1 values used.
  sig = mprage(t1_vals, ttfa, tr, alpha, esp, tf, reps);
  sig = reshape(squeeze(sig), numel(t1_vals) * numel(alpha), tf).';

  if (bias == true)
    vec = ones(tf, 1);
    vec = vec/norm(vec(:));  
    sig = sig - vec * vec' * sig;
    vec = reshape(vec, 1, 1, 1, 1, 1, tf, 1);
    tk  = tk - 1;
  end

  [u, sv, v] = svd(sig, 'econ');
  sv = diag(sv);
  u = u(:, 1:tk);
  u = reshape(u, 1, 1, 1, 1, 1, tf, tk);

  if (bias == true)
    u = cat(7, vec, u);
  end
end

function sig = mprage(t1_vals, ttfa, tr, alpha, esp, tf, reps)
  % function sig = mprage(t1_vals, ttfa, tr, alpha, esp, tf, reps)
  %
  % Returns MPRAGE signal evolution based on input parameters.
  % Assumes inversion.
  %
  % Inputs:
  %   t1_vals  Array of T1 values in milliseconds.
  %   ttfa     Time from 180* to first alpha in milliseconds.
  %   alpha    Small flip angle in degrees.
  %   esp      Echo spacing, or time between adjacent alphas 
  %            in milliseconds.
  %   tf       Turbo factor, or the number of echos.
  %
  % Output:
  %   sig      Signal over time. Dimensions: [size(t1_vals), 1, 1, tf].
  for k=1:numel(size(t1_vals))
    if (k > 3)
      assert(size(t1_vals, k) == 1);
    end
  end
  [sx, sy, sz] = size(t1_vals);
  sig = zeros(sx * sy * sz, numel(alpha), tf);
  for tdx=1:numel(t1_vals)
    for adx=1:numel(alpha)
      t1 = t1_vals(tdx);
      fa = alpha(adx);
      sig(tdx, adx, :) = evol_diffeq(t1, ttfa, tr, fa, esp, tf, reps);
    end
  end
  sig = reshape(sig, sx, sy, sz, numel(alpha), tf);
end

function Mxy = evol_diffeq(t1, ttfa, tr, alpha, esp, tf, reps)
  % function Mxy = evol(t1, ttfa, tr, alpha, esp, tf)
  %
  % Simulates MPRAGE evolution. Helper function to mprage.
  % Please see mprage for input description.

  % Converting to seconds and radians.
  t1    = t1    * 1e-3;
  ttfa  = ttfa  * 1e-3;
  tr    = tr    * 1e-3;
  alpha = alpha * pi/180;
  esp   = esp   * 1e-3;

  M0   = 1;
  Mz   = -M0;
  f    = @(t, Mz) (M0 - Mz) * (1 - exp(-t/t1)) + Mz;
  Mxy  = zeros(tf, 1);

  for p =1:reps
    time = ttfa;
    for q = 1:tf
      Mz     = f(time, Mz);
      Mxy(q) = sin(alpha) * Mz;
      Mz     = cos(alpha) * Mz;
      time   = esp;
    end
    Mz = f(time, Mz);

    time = tr - (ttfa + (tf * esp));
    Mz   = -1 * f(time, Mz);
  end
end
