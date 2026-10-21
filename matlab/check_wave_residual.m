function [rel, R] = check_wave_residual(eta, phi, T)
% CHECK_WAVE_RESIDUAL  How well does a solution satisfy the water wave
% equations of Wilkening's formulation, WITHOUT evaluating the
% Dirichlet--Neumann operator?
%
%   [rel, R] = check_wave_residual(eta, phi, T)
%
%   eta, phi : (Mt x NX) surface elevation and surface potential on the
%              uniform periodic x-grid (NX points on (-pi, pi]) and the
%              uniform time grid t = linspace(0, T/4, Mt). Works equally on
%              Wilkening's stored data (mirrored to 512 points) and on FNO
%              predictions (Y_pred(k, :, 1, :) and Y_pred(k, :, 2, :)).
%   T        : the period (dataset units, g = 1).
%
%   rel      : scalar relative residual (RMS of R over space-time, divided
%              by the RMS size of the equation's terms).
%   R        : (Mt x NX) pointwise residual field.
%
% METHOD. The governing system is
%     eta_t = G(eta) phi                                   (kinematic)
%     phi_t = -g eta - phi_x^2/2
%             + (G(eta)phi + eta_x phi_x)^2 / (2(1+eta_x^2))   (dynamic)
% Evaluating G requires solving Laplace's equation; but the kinematic
% equation itself supplies G(eta)phi = eta_t, so substituting it into the
% dynamic equation leaves a residual computable from the fields alone:
%     R = phi_t + g eta + phi_x^2/2
%         - (eta_t + eta_x phi_x)^2 / (2(1+eta_x^2)).
% x-derivatives are spectral (uniform periodic grid); t-derivatives are
% 4th-order finite differences (2nd order at the two end slices).
%
% WHAT IT DOES AND DOES NOT TEST. R == 0 verifies the Bernoulli (dynamic)
% equation under the kinematic identification of the normal velocity; the
% one property it cannot see is harmonicity itself (that eta_t truly equals
% the Laplace-extension normal derivative of phi). Checking that last link
% requires an independent evaluation of G (Craig--Sulem series, a conformal
% map, or the boundary-integral evaluation from the stored mu).
%
% Calibration on this dataset: Wilkening's stored solutions give
% rel ~ 2-3e-6 at the full 201 timeslices and ~1.5e-5 at 101 slices (the
% floor is set by the time differencing); trained FNO predictions land
% around 6-8e-3. For exploratory predictions with unknown period, minimize
% rel over T to obtain a PREDICTED period.

[Mt, NX] = size(eta);
g  = 1.0;
dt = (T/4) / (Mt - 1);

% spectral d/dx on the uniform periodic grid
kw = [0:NX/2, -NX/2+1:-1] * (2*pi / (2*pi));      % wavenumbers (domain 2*pi)
Dx = @(F) real(ifft(1i .* kw .* fft(F, [], 2), [], 2));

% 4th-order d/dt (interior), 2nd-order one-sided at the ends
Dt = @(F) time_derivative(F, dt);

eta_x = Dx(eta);  phi_x = Dx(phi);
eta_t = Dt(eta);  phi_t = Dt(phi);

R = phi_t + g*eta + 0.5*phi_x.^2 ...
    - (eta_t + eta_x .* phi_x).^2 ./ (2*(1 + eta_x.^2));

scale = sqrt(mean(phi_t(:).^2) + mean(eta(:).^2));
rel   = sqrt(mean(R(:).^2)) / scale;
end

function D = time_derivative(F, dt)
    D = zeros(size(F));
    D(3:end-2, :) = (F(1:end-4, :) - 8*F(2:end-3, :) ...
                   + 8*F(4:end-1, :) - F(5:end, :)) / (12*dt);
    D(1, :)   = (-3*F(1, :) + 4*F(2, :) - F(3, :)) / (2*dt);
    D(2, :)   = (F(3, :) - F(1, :)) / (2*dt);
    D(end-1, :) = (F(end, :) - F(end-2, :)) / (2*dt);
    D(end, :) = (3*F(end, :) - 4*F(end-1, :) + F(end-2, :)) / (2*dt);
end
