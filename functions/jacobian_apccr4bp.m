function A = jacobian_apccr4bp(tau, x_vec, params)
%==========================================================================
%
% Computes Jacobian matrix A = df/dx for autonomous P
% CCR4BP Hamiltonian system.
% Used in STM propagation: dot(Phi) = A*Phi.
%
% Author: G. Montseny
% Date: June , 2026
%
%
% INPUT:               Description                                   Units
%
%  t         -   time                                               -
%  x_vec     -   state vector (6x1)                                 -
%  params    -   parameter struct                                   -
%
% OUTPUT:
%
%  A         -   Jacobian matrix (6 x 6)                              -
%==========================================================================
    
    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    t = x_vec(3);

    switch params.frame
        case 'synodic_m1m2'
            mu = params.mu;
            mu_3 = params.mu_pert;
            r_13 = params.r_pert;
            Omega_3 = params.Omega_pert;
            theta_3_0 = params.theta_pert_0;
            theta_3_t = theta_pert(t, params);
            
            x_1 = -mu; y_1 = 0; r_1_vec = [x_1; y_1];
            x_2 = 1-mu; y_2 = 0; r_2_vec = [x_2; y_2];

            x_3 = -mu + r_13*cos(theta_3_t);
            y_3 = r_13*sin(theta_3_t);
            r_3_vec = [x_3; y_3];

        case 'synodic_m1m3'
        otherwise
            error('Invalid frame')
    end


    % Calculate important quantities
    r_vec = [x; y];
    r_34 = norm(r_vec-r_3_vec);
    
    % Calculate some elements
    H_xx = - (1-mu)*drm1_d2('xx', x, y, - mu, 0)...
        - mu*drm1_d2('xx', x, y, 1- mu, 0)...
        -mu_3*drm1_d2('xx', x, y, x_3, y_3);
    
    H_yy = - (1-mu)*drm1_d2('yy', x, y, - mu, 0)...
        - mu*drm1_d2('yy', x, y, 1- mu, 0)...
        -mu_3*drm1_d2('yy', x, y, x_3, y_3);

    H_xy = - (1-mu)*drm1_d2('xy', x, y, - mu, 0)...
        - mu*drm1_d2('xy', x, y, 1- mu, 0)...
        -mu_3*drm1_d2('xy', x, y, x_3, y_3);

    H_xt = mu_3 * (Omega_3 - 1) *  ( ...
        r_13 * sin(theta_3_t) / r_34^3 ...
        - 3 * r_13 * (x + mu - r_13 * cos(theta_3_t)) * (1 / r_34^5) * ( ...
        sin(theta_3_t) * (x + mu - r_13 * cos(theta_3_t)) - ...
        cos(theta_3_t) * (y - r_13*sin(theta_3_t))) ...
        - sin(theta_3_t) / r_13^2);

    H_yt = mu_3 * (Omega_3 - 1) *  ( ...
        - r_13 * cos(theta_3_t) / r_34^3 ...
        - 3 * r_13 * (y - r_13*sin(theta_3_t)) * (1 / r_34^5) * ( ...
        sin(theta_3_t) * (x + mu - r_13 * cos(theta_3_t)) - ...
         cos(theta_3_t) * (y - r_13*sin(theta_3_t))) ...
        + cos(theta_3_t) / r_13^2);

    H_tt = - 3 * mu_3 * r_13^2 * (Omega_3 - 1)^2 * (1 / r_34^5) * ( ...
        sin(theta_3_t) * (x + mu - r_13*cos(theta_3_t)) ...
        - cos(theta_3_t) * (y - r_13 * sin(theta_3_t)))^2 ...
        + mu_3 * r_13 * (Omega_3-1)^2 * (1/r_34^3) * ( ...
        cos(theta_3_t) * (x + mu) + y * sin(theta_3_t)) ...
        - mu_3 * (Omega_3 -1)^2 * (1 / r_13^2) * ( ...
            x * cos(theta_3_t) + y * sin(theta_3_t));


    % Complete matrix A
    A = [0, 1, 0, 1, 0, 0;
        -1, 0, 0, 0, 1, 0;
        0, 0, 0, 0, 0, 0;
        -H_xx, -H_xy, -H_xt, 0, 1, 0;
        -H_xy, -H_yy, -H_yt,  -1, 0, 0;
        -H_xt, -H_yt, -H_tt, 0, 0, 0];


end