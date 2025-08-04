function solve_three_gamete_genotype_ode()
    % Define parameter values
    mu_val = 1e-4;
    nu_val = 1e-4;
    h1_val = 1;
    h2_val = 1;
    h3_val = 1;
    s_val = 8e-4;
    
    % Initial conditions for original variables
    g0_init = 0.04;
    g1_init = 0.32;
    g2_init = 1 - g0_init - g1_init;
    
    % Initial conditions for log-transformed variables
    L0_init = log(g0_init);
    L1_init = log(g1_init);
    L2_init = log(g2_init);
    
    % Determine the rate limiting parameter
    min_param = min([mu_val, nu_val, s_val]);
    
    % Time span
    tspan = [0 (1/min_param)*10^2];
    
    % Solve the ODE system
    [t, y] = ode15s(@(t, y) odefunc(t, y, s_val, mu_val, nu_val, h1_val, h2_val, h3_val), ...
                    tspan, [L0_init, L1_init, L2_init]);
    
    % Convert back to original variables
    g0 = exp(y(:,1));
    g1 = exp(y(:,2));
    g2 = exp(y(:,3));
    
    % Calculate allele frequencies
    p = g0 + 0.5*g1;  % frequency of ancestral allele
    q = 0.5*g1 + g2;  % frequency of derived allele
    
    % Calculate genotype frequencies for display
    G0 = g0.^2;
    G1 = 2*g0.*g1;
    G2 = 2*g0.*g2 + g1.^2;
    G3 = 2*g1.*g2;
    G4 = g2.^2;
    
    % Plot allele frequencies
    figure;
    subplot(3,1,1);
    hold on;
    plot(t, p, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'p (ancestral)');
    plot(t, q, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'q (derived)');
    xlabel('Time', 'FontSize', 16);
    ylabel('Allele Frequencies', 'FontSize', 16);
    title('Evolution of Allele Frequencies (Autotetraploid Model)', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
    ylim([0 1]);
    
    % Plot gamete frequencies
    subplot(3,1,2);
    hold on;
    plot(t, g0, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'g0');
    plot(t, g1, 'Color', '#A23B72', 'LineWidth', 2, 'DisplayName', 'g1');
    plot(t, g2, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'g2');
    xlabel('Time', 'FontSize', 16);
    ylabel('Gamete Frequencies', 'FontSize', 16);
    title('Evolution of Gamete Frequencies', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
    
    % Plot genotype frequencies
    subplot(3,1,3);
    hold on;
    plot(t, G0, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'G0');
    plot(t, G1, 'Color', '#5B9BD5', 'LineWidth', 2, 'DisplayName', 'G1');
    plot(t, G2, 'Color', '#A23B72', 'LineWidth', 2, 'DisplayName', 'G2');
    plot(t, G3, 'Color', '#F18F01', 'LineWidth', 2, 'DisplayName', 'G3');
    plot(t, G4, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'G4');
    xlabel('Time', 'FontSize', 16);
    ylabel('Genotype Frequencies', 'FontSize', 16);
    title('Evolution of Genotype Frequencies', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
end

function dydt = odefunc(~, y, s, mu, nu, h1, h2, h3)
    % Extract log-transformed variables
    L0 = y(1);
    L1 = y(2);
    L2 = y(3);
    
    % Convert to original gamete variables
    g0 = exp(L0);
    g1 = exp(L1);
    g2 = exp(L2);
    
    % Calculate genotype frequencies from gamete frequencies (Hardy-Weinberg)
    G0 = g0^2;                    % frequency of genotype with 0 derived alleles
    G1 = 2*g0*g1;                 % frequency of genotype with 1 derived allele
    G2 = 2*g0*g2 + g1^2;          % frequency of genotype with 2 derived alleles
    G3 = 2*g1*g2;                 % frequency of genotype with 3 derived alleles
    G4 = g2^2;                    % frequency of genotype with 4 derived alleles
    
    % Mean fitness calculation
    wbar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4);
    
    % Relative fitness calculations
    w0 = 1/wbar;
    w1 = (1 - s*h1)/wbar;
    w2 = (1 - s*h2)/wbar;
    w3 = (1 - s*h3)/wbar;
    w4 = (1 - s)/wbar;
    
    % Selection terms
    sel_g0 = G0*w0 + (1/2)*G1*w1 + (1/6)*G2*w2;
    sel_g1 = (1/2)*G1*w1 + (2/3)*G2*w2 + (1/2)*G3*w3;
    sel_g2 = (1/6)*G2*w2 + (1/2)*G3*w3 + G4*w4;
    
    % Mutation terms
    mut_g0 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2) - g0;
    mut_g1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu) + 2*sel_g2*(1-nu)*nu - g1;
    mut_g2 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2) - g2;
    
    % Note: The constraint g0 + g1 + g2 = 1 should be maintained automatically
    % by the structure of the equations, but we could enforce it by setting:
    % mut_g2 = -mut_g0 - mut_g1;
    
    % Log-transformed ODEs
    dL0dt = (1 / g0) * mut_g0;
    dL1dt = (1 / g1) * mut_g1;
    dL2dt = (1 / g2) * mut_g2;
    
    dydt = [dL0dt; dL1dt; dL2dt];
end

% Run the simulation
solve_three_gamete_genotype_ode