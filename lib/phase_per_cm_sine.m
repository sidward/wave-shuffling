function [ppcm, G0] = phase_per_cm_sine(wx, adc, num_cyc, gmax, smax, dt, sc)
  % [ppcm, G0] = phase_per_cm_sine(wx, adc, num_cyc, smax, smax, dt)
  %
  % Generates to phase per cm for a sine gradient wave.
  %
  % Inputs:
  %   wx            Readout dimension after taking into account aliasing due to wave.
  %   adc           Readout time in ms.
  %   num_cyc       Number of cycles.
  %   gmax          Max gradient amplitude in Gauss per cm.
  %   smax          Max gradient slew rate in Gauss per cm per second.
  %   dt            Gradient delay. (Optional.)
  %   sc            Gradient scale. (Optional.)
  %
  % Output:
  %   ppcm          Phase per cm that can be used by wave_psf. 
  %   G0            Actually gradient amplitude used.
  %
  % Original written by Kawin Setsompop.
  % Modified as needed.

  if (nargin < 6)
    dt = 0;
    sc = 1;
  end

  ADCduration = adc * 1e3;

  wavepoints = floor(ADCduration/10);
  T_wavepoints = 10^-5;
  TimePerSine = (wavepoints) * T_wavepoints/num_cyc; 
  w = 2 * pi / TimePerSine;

  if w * smax >= gmax
    G0 = gmax;   % Gradient limited.
  else
    G0 = smax/w; % Slew limited.
  end

  GradTimePoints = ([1:wavepoints] * T_wavepoints);
  Gwave = sc * G0 * sin(w * ((GradTimePoints - T_wavepoints) - dt));
  Gwave = [0, Gwave]; GradTimePoints = [0, GradTimePoints];
  GwaveforPhase = ([0 Gwave] + [Gwave 0])/2; GwaveforPhase = GwaveforPhase(1:end-1);
  phasePerCm = 2 * pi * cumsum(GwaveforPhase) * 4257.56 * T_wavepoints;

  PrePhase   = -(G0/w) * 2 * pi * 4257.56;
  phasePerCm = phasePerCm + PrePhase;

  % Resampling to get on kx sampling rate.
  T_KxSampling = (ADCduration * 1E-6)/wx;
  KxTimePoints = [1:wx] * T_KxSampling;
  ppcm = interp1(GradTimePoints, phasePerCm, KxTimePoints, 'spline');
end
