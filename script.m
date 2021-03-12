%% MPRAGE WSHFL
%
clear all; close all; clc;

%% Paths ----------------------------------------------------------------------------------------------------
addpath(genpath('lib'));

%% Parameters -----------------------------------------------------------------------------------------------

% Data dimensions.
SX   = 256;                 % Readout dimension.
SY   = 256;                 % Phase encode 1.
SZ   = 256;                 % Phase encode 2.
sy   = 512;                 % Phase encode 1 after zero-pad for linear convolution.
sz   = 512;                 % Phase encode 2 after zero-pad for linear convolution.
nc   = 16;                  % Number of channels after coil compression.
tf   = 256;                 % Echo train length.
tk   = 2;                   % Shuffling temporal basis rank.
OS   = 6;                   % Readout over-sampling factor for wave.
nadc = 8192;                % Number of ADCs acquired.
wx   = SX * OS;             % Readout dimension after zero-pad for wave-PSF.

% Temporal parameters.
esp   = 8.4;                % Echo spacing in ms.
TR    = 2500;               % Repition time in ms.
TI    = 1100;               % Inversion time in ms.
ttfa  = TI - tf/2 * esp;    % Time between 180* and first alpha.
alpha = 9;                  % Flip angle of excitation in degrees.

% Wave parameters.
dy  = 0.1 * (SY/sy);        % Resolution in cm.
dz  = 0.1 * (SZ/sz);        % Resolution in cm.
adc = 5.0688;               % ADC duration in ms.
gm  = 2.7;                  % Max gradient amplitude in Gauss/cm.
sm  = 18700.00;             % Gradient slew rate in Gauss/cm/s.
cyc = 5;                    % Number of wave-cycles.
ofY = -3.5;                 % Offset in Y in cm.
ofZ = 0.05;                 % Offset in Z in cm.

% ShflPSF parameters.
sl  = 100;                  % Sigmoid lambda.
tol = 1E-8;                 % Tolerance value.

%% Download data --------------------------------------------------------------------------------------------
if (~isfile('data/tbl.cfl'))
  fprintf('Downloading data. ---------------------------------------------------------------------------\n');
  fprintf('If certificate errors, edit `download.sh` to `wget --no-check-certificate ...`.\n');
  system('cd data/ && bash download.sh && cd ../');
end

%% Create espirit maps and coil compression matrix ----------------------------------------------------------
if (~isfile('data/whiten_cc_mat.cfl'))
  fprintf('Preparing reference and whitening matrix. ---------------------------------------------------\n');

  % Load data.
  rref = readcfl('data/reference_scan');
  nmat = readcfl('data/noise_scan');

  % Account for over-sampling in reference scan.
  rref = F_inv(rref, [1, 2, 3]);
  rref = rref(floor((size(rref, 1) * (1 - 1/OS))/2) + [1:size(rref, 1)/OS], :, :, :);
  [rx, ry, rz, rc] = size(rref);

  % SVD to get coil compression matrix.
  rref = reshape(rref, rx * ry * rz, rc).'; 
  [u, s, ~] = svd(rref, 'econ');
  u = u(:, 1:nc);
  s = diag(s);

  % Coil compressing noise.
  nmat = u' * nmat;

  % Estimating whitening matrix.
  covm = (nmat * nmat')/(size(nmat, 2) - 1);
  whmt = inv(chol(covm, 'lower'));

  % Joint coil compression and whitening matrix.
  whtcc = whmt * u';

  % Whiten and coil compress reference data.
  rref = whtcc * rref;
  rref = rref.';
  rref = reshape(rref, rx, ry, rz, nc);
  rref = F_fwd(rref, [1, 2, 3]);

  % Resizing data for BART espirit.
  ref = zeros(SX, sy, sz, nc);
  ref(floor((SX - size(rref, 1))/2) + [1:size(rref, 1)], ...
      floor((sy - size(rref, 2))/2) + [1:size(rref, 2)], ...
      floor((sz - size(rref, 3))/2) + [1:size(rref, 3)], ...
      :) = rref;

  % Saving data.
  writecfl('data/reference',     ref);
  writecfl('data/whiten_cc_mat', whtcc);
end

%% Load data ------------------------------------------------------------------------------------------------
if (~isfile('data/tblcc.cfl'))
  fprintf('Preparing ADC table. ------------------------------------------------------------------------\n');
  whtcc = readcfl('data/whiten_cc_mat');
  tbl   = readcfl('data/tbl');

  tbl  = permute(tbl, [2, 1, 3]);
  rc   = size(tbl, 1);
  nacq = size(tbl, 3);

  tbl = reshape(tbl, rc, wx * nacq);
  tbl = reshape(whtcc * tbl, nc, wx, nacq);
  tbl = permute(tbl, [2, 1, 3]);

  writecfl('data/tblcc', tbl);
end

%% Temporal basis -------------------------------------------------------------------------------------------
if (~isfile('data/phi.cfl'))
  fprintf('Preparing temporal basis. -------------------------------------------------------------------\n');
  alpha = [0.7:0.1:1.3] * alpha;
  [phi, sig, sv, t1v] = mprage_basis(2 * tk, [50:10:5000], ttfa, TR, alpha, esp, tf, 100, false);
  writecfl('data/initphi', reshape(phi(:, :, :, :, :, :, 1), 1, 1, 1, 1, 1, tf, 1));
  writecfl('data/phi',     phi(:, :, :, :, :, :, 1:tk));
end

%% Calibrate maps -------------------------------------------------------------------------------------------
if (~isfile('data/maps.cfl'))
  fprintf('ESPIRiT. ------------------------------------------------------------------------------------\n');
  system(sprintf('bart ecalib -v 1 -m 1 -a -P %s %s', ...
                    'data/reference',                 ...
                    'data/maps'));
end

%% Generate PSF ---------------------------------------------------------------------------------------------
if (~isfile('data/initpsf.cfl'))
  fprintf('Preparing initial wave-psf. -----------------------------------------------------------------\n');

  yaxis = dy * ((0:(sy-1)) - floor(sy/2)) - ofY;
  zaxis = dz * ((0:(sz-1)) - floor(sz/2)) - ofZ;

  pY  = phase_per_cm_sine  (wx, adc, cyc, gm, sm);
  pZ  = phase_per_cm_cosine(wx, adc, cyc, gm, sm);

  fl = @(x) x(:);

  cY = fft(pY);
  cY(1 + (end/2):end) = 0;
  [~, lY] = sort(abs(fl(cY)), 'descend');
  lY = lY(1);
  cY = [fl(real(cY(lY))); fl(imag(cY(lY)))];

  cZ = fft(pZ);
  cZ(1 + (end/2):end) = 0;
  [~, lZ] = sort(abs(fl(cZ)), 'descend');
  lZ = lZ(1);
  cZ = [fl(real(cZ(lZ))); fl(imag(cZ(lZ)))];

  c = [cY(:); cZ(:)];
  l = [lY(:); lZ(:)];
  psf = wave_coeffs_to_psf(c, l, wx, yaxis, zaxis);

  writecfl('data/initpsf',    psf);
  writecfl('data/initcoeffs', c);
  writecfl('data/initlocs',   l);
end

%% Initial rec ----------------------------------------------------------------------------------------------
if (~isfile('data/initrec.cfl'))
  fprintf('Initial reconstruction. ---------------------------------------------------------------------\n');
  system(sprintf('bart wshfl -H -i 50 -R W:7:0:0.005 %s %s %s %s %s %s', ...
    'data/maps',      ...
    'data/initpsf',   ...
    'data/initphi',   ...
    'data/rdr',       ...
    'data/tblcc',     ...
    'data/initrec'));
end

%% Calibrate PSF --------------------------------------------------------------------------------------------
if (~isfile('data/wave.cfl'))
  fprintf('ShflPSF. ------------------------------------------------------------------------------------\n');
  w = readcfl('data/initpsf');
  x = squeeze(readcfl('data/initrec'));
  c = squeeze(readcfl('data/initcoeffs'));
  l = squeeze(readcfl('data/initlocs'));

  y_axis = dy * ((0:(sy-1)) - floor(sy/2)) - ofY;
  z_axis = dz * ((0:(sz-1)) - floor(sz/2)) - ofZ;

  [w, c, shift, locs, res] = shflpsf(x, w, y_axis, z_axis, c, l, 1, sl, tol);

  writecfl('data/wave', w);
end

%% Final reconstruction -------------------------------------------------------------------------------------
if (~isfile('coeffs.cfl'))
  fprintf('Final reconstruction. -----------------------------------------------------------------------\n');
  system(sprintf('bart wshfl -H -i 200 -R W:7:0:0.002 %s %s %s %s %s %s', ...
    'data/maps',   ...
    'data/wave',   ...
    'data/phi',    ...
    'data/rdr',    ...
    'data/tblcc',  ...
    'coeffs'));
end

%% Time series ----------------------------------------------------------------------------------------------
if (~isfile('time_series.cfl'))
  fprintf('Preparing time-series. ----------------------------------------------------------------------\n');
  fprintf('NOTE: For memory reasons, only saving a single slice.\n');
  phi = readcfl('data/phi');
  img = readcfl('coeffs');
  img = img(104, 1:2:end, 1:2:end, :, :, :, :);
  ts  = abs(sum(img .* phi, 7));
  writecfl('time_series', ts);
end

%% ----------------------------------------------------------------------------------------------------------
