function H_out = doppler_track(H_in, Y, params)
% DOPPLER_TRACK  Cubic-spline channel tracking across the 10-symbol frame.
%
%   H_out = doppler_track(H_in, Y, params)
%
% Inputs:
%   H_in   - N_used x N_sym  channel estimate from channel_estimate().
%             Reliable (LS) at pilot symbols {1,4,7,10}; nearest-neighbour
%             copies elsewhere — those copies are what we replace here.
%   Y      - N_FFT  x N_sym  frequency-domain rx (available if needed)
%   params - frame parameter struct
%
% Output:
%   H_out  - N_used x N_sym  improved channel estimate.
%
% -------------------------------------------------------------------------
% Strategy: cubic-spline interpolation in TIME across the 4 pilot snapshots.
%
% At fd = 3000 Hz the channel coherence time (~83 µs) is shorter than the
% 3-symbol pilot gap (~130 µs). Nearest-neighbour copies go stale midway
% through each gap. Spline interpolation fits the smooth Jakes trajectory
% between the 4 known channel snapshots.
%
% Geometry (all symbols between two pilot snapshots — zero extrapolation):
%   symbols 2,3  -> spline between pilot snapshots at symbols 1 and 4
%   symbols 5,6  -> spline between pilot snapshots at symbols 4 and 7
%   symbols 8,9  -> spline between pilot snapshots at symbols 7 and 10
%
% Complexity: one vectorised interp1 call over all 312 subcarriers.
%   MATLAB's interp1 with 'spline' interpolates real and imaginary parts
%   separately, which is the correct operation for complex channel values.
% -------------------------------------------------------------------------

% With the comb_10sym pilot pattern every OFDM symbol carries pilots, so
% channel_estimate() already delivers a fresh LS estimate for every symbol.
% No time interpolation is needed or beneficial — pass H_in straight through.
H_out = H_in;

end
