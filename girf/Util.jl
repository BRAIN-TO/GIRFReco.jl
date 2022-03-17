## Authors:   Johanna Vannesjo (johanna.vannesjo@gmail.com),
#             Alexander Jaffray (alexander.jaffray@gmail.com),
#             Lars Kasper,
#             Tim Zhe Wu
#
# Copyright (C) 2014 IBT, University of Zurich and ETH Zurich,
#               2016 FMRIB centre, University of Oxford,
#               2021 BRAIN-TO, University of Toronto
#
# This file is part of a code package for GIRF computation and application.
# The package is available under a BSD 3-clause license. Further info see:
# https://github.com/MRI-gradient/girf

## Port of Johanna Vannesjo's time2freq method
#  https://github.com/MRI-gradient/GIRF/blob/main/utils/time2freq.m

function time2freq(t)

    nrs = length(t)
    dt = t[2] .- t[1]
    f_max = 1 ./dt
    df = f_max ./ nrs
    f = ((0:nrs-1) .- floor(nrs/2.0)).*df

    return f, df, f_max

end

## BW_filter applies a frequency domain filter to the data
#
# The Function BW_filter applies a frequency-domain filter to the data.
# Typical input-data are k-evolutions that can be filtered given the finite
# bandwith of the gradient chain.
# The filter BW corresponds to the FWHM. Three different filters are
# disposable. The raised cosine filter is set as default. However, the filter BW makes much more % of a difference than
# the filter-type.
#
# USE
# [k_data] = BW_filter(k_data, t, BW, filterType, beta)
#
# IN
#   k_data      [nr_samples nr_orders nDim3]
#   t           [nr_samples 1] in [sec]
#   BW          [scalar] in [Hz] bandwidth at FWHM of the filter
#   filterType  [char], optional: 'ga', 'bl', 'bh' or 'rc': the filter type (gaussian, blackman, blackman-harris or raised cosine)
#   beta_rc     [scalar], optional: roll-off factor between 0 and 1 for the raised cosine filter (0 giving a box-car function, and 1 a cosine without plateau)
#
# OUT
#   k_data      [nr_samples x nr_orders]; filtered data
#
# TODO

function BW_filter(k_data, t, BW, filterType, beta)

  @error "NOT IMPLEMENTED YET"

end

## Applies selected bandwidth window to input array
#
# The Function BW_window applies a frequency-domain windowing to the data.
# The filter BW corresponds to the FWHM. Available filters are Gaussian,
# Blackman, Blackmanharris and Raised cosine.
#
# IN:
#   array       [nr_samples x nr_orders] data to be windowed
#   f           [nr_samples x 1] in [Hz] frequency vector of the array
#   BW          [scalar] in [Hz] bandwidth at FWHM of the filter
#   filterType  {'ga' 'bl' 'bh' 'rc'} the filter type (gaussian, blackman, blackman-harris or raised cosine)
#   beta        [scalar] roll-off factor for raised cosine filter
#
# OUT:
#   array       [nr_samples x nr_orders] windowed data
#   filter      [nr_samples x 1] calculated filter
#
# TODO

function BW_window(array, f, BW, filterType, beta)

  @error "NOT IMPLEMENTED YET"

end

## Function to compute unwrapped phase of complex vector, with zero phase at zero
# frequency (assuming fftshift)
#
# USE
# phase = CenteredPhase(array)
#
# IN
# array     [n_samples x n_columns] complex-valued array
#
# OUT
# phase     [n_samples x n_columns] phase of array
#
# TODO

function CenteredPhase(array)

  @error "NOT IMPLEMENTED YET"

end

## Function to compute time vector, with zero at center (assuming fftshift)
#
# USE
# t = CenteredTime(dt, nrs)
#
# IN
# dt    sampling interval
# nrs   number of samples
#
# OUT
# t     time vector, centered at zero (assuming fftshift)
#
# TODO

function CenteredTime(dt, nrs)

  @error "NOT IMPLEMENTED YET"

end

## Function to compute set of blip inputs
#
# TODO

function ComputeInputs()

  @error "NOT IMPLEMENTED YET"

end

## Function creates a raised cosine window for frequency domain filtering
#
# USE
# [filter] = raised_cosine(f, T, beta)
#
# IN
#   f       [nr_samples x 1]; frequency vector
#   T       [scalar]; 1/BW at FWHM of the window
#   beta    [scalar]; roll-off factor between 0 and 1 for the raised-cosine
#                     window (0 giving a box-car function, and 1 a cosine without plateau)
#
# OUT
#   filter  [nr_samples x 1]; raised cosine window
#
# See also BW_filter
#
# TODO

function raised_cosine(f, T, beta)

  @error "NOT IMPLEMENTED YET"

end

## Function to create frequency-swept pulses from parameters
#
# USE
# sweep = sweeps(t, T_acq, f1, f2, phi0, A, t_start, type)
#
# IN
# t         [n_samples x 1] time vector [s]
# T_acq     length of sweep [s]
# f1        starting frequency
# f2        end frequency
# phi0      starting phase
# A         pulse amplitude
# t_start   pulse start time
# type      which frequency-modulation of pulse?
#
# OUT
# sweep     [n_samples x 1] frequency-swept pulse
# f_t       [n_samples x 1] frequency-modulation as function of time
#
# TODO

function sweeps(t, T_acq, f1, f2, phi0, A, t_start, type)

  @error "NOT IMPLEMENTED YET"

end

## Function to create trapezoidal pulse shapes from parameters
#
# USE
# grads = trapezoid(t,ons,amp,dur,plateau,dur2)
#
# IN
# t         [n_samples x 1] time vector [s]
# ons       [1 x n_pulses] onset of pulse [s]
# amp       [1 x n_pulses] amplitude of pulse
# dur       [1 x n_pulses] duration of first slope
# plateau   [1 x n_puless] duration of pulse plateau (opt: varargin{1})
# dur2      [1 x n_pulses] duration of second slope (opt: varargin{2})
#
# OUT
# grads     [n_samples x n_pulses] trapezoidal gradient pulse
#
# TODO

function trapezoid(t,ons,amp,dur,plateau,dur2)

  @error "NOT IMPLEMENTED YET"

end

## Performs frequency-dependent smoothing of GIRF/SIRF
#
# USE
# [SIRF] = VariableSmoothing(SIRF, f, f2)
#
# IN
#   SIRF    [nr_samples nr_k] SIRF matrix to be smoothed
#   f       [nr_samples 1] frequency vector
#   fMax    [Hz] maximum frequency to be smoothed & cut
#   BWrange [BWmin BWmax] [Hz] minimal and maximal bandwidth of smoothing kernel (optional)
#   vsBW    [Hz] frequency region of narrower smoothing kernel (optional)
#
# OUT
#   SIRF    [nr_samples nr_k] smoothed SIRF
#   fs      [nr_samples 1] cut frequency vector
#
# TODO

function VariableSmoothing(SIRF, f, f2)

  @error "NOT IMPLEMENTED YET"

end


# ## Currently only works for spiral trajectories but will have to be extended to cartesian!
# function parseTrajectoryGradients(a::AcquisitionData)

#   for l = 1:length(a.traj)
      
#       nProfiles = a.traj[l].numProfiles
#       nSamples = a.traj[l].numSamplingPerProfile
#       nodes = a.traj[l].nodes

#       for profile = 1:nProfiles
          
#           ilExtractor = nSamples*(profile-1) .+ (1:nSamples)
#           ilNodes = nodes[:,ilExtractor]

#           figure()
#           ilGrads = nodes_to_gradients(ilNodes)

#           plot(ilGrads')

#       end
  
#   end
      
# end

## Get gradients from the trajectory
function nodes_to_gradients(nodes::Matrix; gamma=42577478, dwellTime=2e-6, FOV=[220,220,1],reconSize=[200,200,1])

  ## Normalized Conversion (norm kspace to grads in T/m) is scalingFactor = reconSize/(gamma*dwellTime*FOV)
  conversionFactor = reconSize./(gamma*dwellTime.*FOV)*1000 # The 1000 factor is the conversion from mm to m

  gradients = diff(hcat([0; 0], nodes), dims = 2)
  gradients = gradients .* conversionFactor[1:2]
  return gradients

end

## Convert gradients to trajectory nodes
function gradients_to_nodes(gradients::Matrix; gamma=42577478, dwellTime=2e-6, FOV=[220,220,1],reconSize=[200,200,1])

  ## Normalized Conversion (grads in T/m to normalized k-space) is scalingFactor = (gamma*dwellTime*FOV)/reconSize
  conversionFactor = ((gamma*dwellTime.*FOV) ./ reconSize) /1000 # The 1000 is the conversion from mm to m

  nodes = cumsum(gradients, dims = 2)
  nodes = nodes .* conversionFactor[1:2]
  return nodes

end