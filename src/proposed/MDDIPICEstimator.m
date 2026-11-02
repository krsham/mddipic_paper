classdef MDDIPICEstimator < AggressiveMDDIPICEstimator
    % MDDIPICESTIMATOR Public interface to joint three-dimensional MDDIPIC.
    %
    % est = MDDIPICEstimator(name,value,...)
    % [tau,nu,theta,g,P] = step(est,yDD,PExact,sigma2DD)
    %
    % The constructor fixes Joint3DRefinement=true. Delay, Doppler and angle are
    % refined together at every shrinking-grid level; all active gains are fitted
    % jointly. Other configuration and diagnostics are inherited from
    % AggressiveMDDIPICEstimator. Units: seconds, Hz, radians, complex amplitude.
    % See capsule_make_estimator for a complete configuration example.
    %
    % SPDX-License-Identifier: GPL-3.0-only
    % Copyright (c) 2026 Koorosh Ahmadi. See COPYRIGHT.txt.
    % Additional MATLAB permission: LICENSES/MATLAB_EXCEPTION.txt.
    % Attribution and paper citation: CITATION.md.

    methods

        function obj = MDDIPICEstimator(varargin)
            names = string(varargin(1:2:end));
            if any(strcmpi(names, "Joint3DRefinement"))
                error('MDDIPICEstimator:CanonicalMode', ...
                      ['Joint3DRefinement is fixed true for canonical ', ...
                       'MDDIPICEstimator. Use LegacySplitMDDIPICEstimator ', ...
                       'for the historical split-refinement behavior.']);
            end
            obj@AggressiveMDDIPICEstimator(varargin{:}, ...
                                           'Joint3DRefinement', true);
        end

    end
end
