function dY_dt_vec = eom_ext_apccr4bp(tau, Y_vec, params)
%==========================================================================
%
% Computes aCCR4BP Hamiltonian equations of motion together with STM dynamics.
% Augmented state is defined as Y = [x; vec(Phi)], where x is 6x1 and Phi is 6x6.
%
% Author: G. Montseny
% Date: June 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  t         -   time                                               -
%  Y_vec     -   augmented state (42x1)                             -
%                 Y_vec = [x(6x1); vec(Phi)(36x1)]
%  params    -   parameter struct                                   -
%
% OUTPUT:
%
%  dY_dt_vec -   time derivative of augmented state (42x1)          -
%==========================================================================

    % Initialization
    Y_vec = Y_vec(:);
   
    % Extract state and STM vector
    x_vec = Y_vec(1:6);
    Phi_vec = Y_vec(7:42);

    % EoM
    dx_dt_vec = eom_apccr4bp(tau, x_vec, params);

    % STM
    Phi_mtx = reshape(Phi_vec, 6, 6);
    A_t = jacobian_apccr4bp(tau, x_vec, params);
    dPhi_dt_mtx = A_t*Phi_mtx;

    % Put vectors back into Y
    dx_dt_vec = dx_dt_vec(:);
    dPhi_dt_vec = dPhi_dt_mtx(:);
    dY_dt_vec = [dx_dt_vec; dPhi_dt_vec];
end