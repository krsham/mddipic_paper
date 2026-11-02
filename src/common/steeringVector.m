function a = steeringVector(theta_rad, Na)
    % STEERINGVECTOR Computes the standard ULA steering vector.
    % Assumes half-wavelength antenna spacing (d = lambda/2).
    %
    % Inputs:
    %   theta_rad - Angle of Arrival/Departure in radians
    %   Na        - Number of antennas
    %
    % Output:
    %   a         - Na x 1 complex steering vector

    % Standard ULA response with d = lambda/2
    % Phase shift between adjacent elements is pi * sin(theta)
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    n_vec = (0:Na - 1).' - Na / 2;
    a = exp(1j * pi * sin(theta_rad) * n_vec);
end
