function solve_dip_ode()
    % Define parameter values
    mu_val = 1e-4;
    nu_val = 1e-4;
    h_val = 0.5;     % dominance coefficient (0.5 = additive, 0 = recessive, 1 = dominant)
    s_val = 8e-4;
    
    % Initial conditions for original variables
    g0_init = 0.1;   % initial frequency of allele 0 (p)
    
    % Initial conditions for log-transformed variable (only need one due to constraint)
    L0_init = log(g0_init);
    
    % Determine the rate limiting parameter
    min_param = min([mu_val, nu_val, s_val]);
    
    % Time span
    tspan = [0 (1/min_param)*10^2];
    
    % Solve the ODE system (only solving for g0, g1 is determined by constraint)
    [t, y] = ode15s(@(t, y) odefunc(t, y, s_val, mu_val, nu_val, h_val), ...
                    tspan, L0_init);
    
    % Convert back to original variables
    g0 = exp(y(:,1));
    g1 = 1 - g0;  % constraint: g0 + g1 = 1
    
    % Calculate genotype frequencies for display
    G0 = g0.^2;           % frequency of genotype 00
    G1 = 2*g0.*g1;        % frequency of genotype 01 (heterozygote)
    G2 = g1.^2;           % frequency of genotype 11
    
    % Plot allele frequencies
    figure;
    subplot(2,1,1);
    hold on;
    plot(t, g0, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'g0 (p)');
    plot(t, g1, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'g1 (q)');
    xlabel('Time', 'FontSize', 16);
    ylabel('Allele Frequencies', 'FontSize', 16);
    title('Evolution of Allele Frequencies (Diploid Model)', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
    ylim([0 1]);
    
    % Plot genotype frequencies
    subplot(2,1,2);
    hold on;
    plot(t, G0, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'G0');
    plot(t, G1, 'Color', '#A23B72', 'LineWidth', 2, 'DisplayName', 'G1');
    plot(t, G2, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'G2');
    xlabel('Time', 'FontSize', 16);
    ylabel('Genotype Frequencies', 'FontSize', 16);
    title('Evolution of Genotype Frequencies', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
    ylim([0 1]);
end

function dydt = odefunc(~, y, s, mu, nu, h)
    % Extract log-transformed variable
    L0 = y(1);
    
    % Convert to original allele variables
    g0 = exp(L0);
    g1 = 1 - g0;  % constraint: g0 + g1 = 1
    
    % Mean fitness calculation
    w_bar = 1 - s*(h*2*g0*g1 + g1^2);
    
    % Relative fitness calculations
    w0 = 1/w_bar;           % fitness of genotype 00
    w1 = (1 - s*h)/w_bar;   % fitness of genotype 01 (heterozygote)
    w2 = (1 - s)/w_bar;     % fitness of genotype 11
    
    % Selection terms (post-selection allele frequencies)
    sel_g0 = w1*g1*g0 + w0*g0^2;
    sel_g1 = w2*g1^2 + w1*g1*g0;
    
    % Mutation terms (change in allele frequency due to selection and mutation)
    mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
    
    % Note: mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1, but we use constraint
    % The constraint g0 + g1 = 1 is maintained automatically since:
    % d/dt(g0 + g1) = mut_g0 + mut_g1 = 0 by construction of the equations
    
    % Log-transformed ODE
    dL0dt = (1 / g0) * mut_g0;
    
    dydt = dL0dt;
end

% Run the simulation
solve_dip_ode