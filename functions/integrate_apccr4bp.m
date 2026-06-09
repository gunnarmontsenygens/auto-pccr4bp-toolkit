function [t_hist, x_vec_hist, Phi_mtx_hist, i_e, t_e, x_e_vec, Phi_e_mtx] = integrate_apccr4bp(tspan, x_0_vec, params, event_fun)
%==========================================================================
%
% Integrates planar autonomous CCR4BP Hamiltonian equations of motion together with the STM.
% Augmented state is defined as X = [x; vec(Phi)].
% All outputs are returned together in the struct traj.
%
% Author: G. Montseny
% Date: June 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  tspan      -   integration interval [t0 tf]                        -
%  x_0_vec    -   initial state (6x1)                                 -
%  params     -   parameter struct                                    -
%  event_fun  -   (optional) event function handle                    -
%
% OUTPUT:
%
%  t_hist        -   time vector                                      -
%  x_vec_hist   -   state history (N x 6)                             -
%  Phi_mtx_hist -   STM history (N x 6 x 6)                           -
%  t_e          -   event times                                       -
%  x_e_vec      -   state at events                                   -
%  i_e          -   event indices                                     -
%  Phi_e_mtx          -   event indices                               -
%==========================================================================
    
    % Check if there's any events
    if nargin < 4
        event_fun = [];
    end
    
    % Extract base ODE options from params 
    ode_options = params.ode.options;
    
    % Implement event in case it's not empy
    if ~isempty(event_fun)
        ode_options = odeset(ode_options, ...
            'Events', @(t,Y) event_fun(t,Y,params));
    end
    
    % Initialization
    x_0_vec = x_0_vec(:);
    Phi_0_mtx = eye(6);
    X_0_vec = [x_0_vec; Phi_0_mtx(:)];
    
    % Integration
    if isempty(event_fun)

        [t_hist, X_vec_hist] = ode113( ...
            @(t,X) eom_ext_apccr4bp(t, X, params), ...
            tspan, X_0_vec, ode_options);

        % Empty event outputs
        t_e = [];
        x_e_vec = [];
        Phi_e_mtx = [];
        i_e = [];

    else

        [t_hist, X_vec_hist, t_e, X_e_vec, i_e] = ode113( ...
            @(t,X) eom_ext_apccr4bp(t, X, params), ...
            tspan, Y_0_vec, ode_options);
        x_e_vec = X_e_vec(1:6);
        Phi_e_mtx = reshape(X_e_vec(7:42),6,6);

    end

    % Recover state history
    x_vec_hist = X_vec_hist(:,1:6);

    % Recover STM history as (N x 6 x 6)
    N = size(X_vec_hist,1);
    Phi_mtx_hist = zeros(N,6,6);
    
    for k = 1:N
        Phi_mtx_hist(k,:,:) = reshape(X_vec_hist(k,7:42),6,6);
    end

end