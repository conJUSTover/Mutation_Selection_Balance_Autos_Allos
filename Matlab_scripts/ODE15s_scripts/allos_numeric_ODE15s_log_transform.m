function solve_allo_ode()
    % Define parameter values
    mu_val = 1e-4;
    nu_val = 1e-4;
    h1_val = 1;
    h2_val = 1;
    h3_val = 1;
    s_val = 8e-4;
    
    % Initial conditions for original variables
    g00_init = 0.6;
    g01_init = 0.15;
    g11_init = 0.1;
    
    % Initial conditions for log-transformed variables
    L00_init = log(g00_init);
    L01_init = log(g01_init);
    L11_init = log(g11_init);
    
    % Determine the rate limiting parameter
    min_param = min([mu_val, nu_val, s_val]);
    
    % Time span
    tspan = [0 (1/min_param)*10^2];
    
    % Solve the ODE system (only solving for 3 variables due to constraint)
    [t, y] = ode15s(@(t, y) odefunc(t, y, s_val, mu_val, nu_val, h1_val, h2_val, h3_val), ...
                    tspan, [L00_init, L01_init, L11_init]);
    
    % Convert back to original variables
    g00 = exp(y(:,1));
    g01 = exp(y(:,2));
    g11 = exp(y(:,3));
    g10 = 1 - g00 - g01 - g11;  % reconstruct g10 from constraint
    
    % Calculate allele frequencies
    pa = g00 + g01;  % frequency of ancestral allele in subgenome a
    pb = g00 + g10;  % frequency of ancestral allele in subgenome b
    qa = g10 + g11;  % frequency of derived allele in subgenome a
    qb = g01 + g11;  % frequency of derived allele in subgenome b
    
    % Calculate genotype frequencies
    G0 = g00.^2;                                      % 0 derived alleles
    G1 = 2*(g00.*g01 + g00.*g10);                     % 1 derived allele
    G2 = g01.^2 + 2*(g00.*g11 + g01.*g10) + g10.^2;   % 2 derived alleles
    G3 = 2*(g01.*g11 + g10.*g11);                     % 3 derived alleles
    G4 = g11.^2;                                      % 4 derived alleles
    
    % Plot results
    figure;
    
    % Top panel: Allele frequencies
    subplot(3,1,1);
    hold on;
    plot(t, pa, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'pa (ancestral a)');
    plot(t, pb, 'Color', '#A23B72', 'LineWidth', 2, 'DisplayName', 'pb (ancestral b)');
    plot(t, qa, 'Color', '#F18F01', 'LineWidth', 2, 'DisplayName', 'qa (derived a)');
    plot(t, qb, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'qb (derived b)');
    xlabel('Time', 'FontSize', 16);
    ylabel('Allele Frequencies', 'FontSize', 16);
    title('Evolution of Allele Frequencies', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
    ylim([0 1]);
    
    % Middle panel: Gamete frequencies
    subplot(3,1,2);
    hold on;
    plot(t, g00, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'g00');
    plot(t, g01, 'Color', '#A23B72', 'LineWidth', 2, 'DisplayName', 'g01');
    plot(t, g10, 'Color', '#F18F01', 'LineWidth', 2, 'DisplayName', 'g10');
    plot(t, g11, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'g11');
    xlabel('Time', 'FontSize', 16);
    ylabel('Gamete Frequencies', 'FontSize', 16);
    title('Evolution of Gamete Frequencies (Four-Gamete Model)', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
    ylim([0 1]);
    
    % Bottom panel: Genotype frequencies
    subplot(3,1,3);
    hold on;
    plot(t, G0, 'Color', '#1E576F', 'LineWidth', 2, 'DisplayName', 'G0 (0 alleles)');
    plot(t, G1, 'Color', '#5B9BD5', 'LineWidth', 2, 'DisplayName', 'G1 (1 allele)');
    plot(t, G2, 'Color', '#A23B72', 'LineWidth', 2, 'DisplayName', 'G2 (2 alleles)');
    plot(t, G3, 'Color', '#F18F01', 'LineWidth', 2, 'DisplayName', 'G3 (3 alleles)');
    plot(t, G4, 'Color', '#C73E1D', 'LineWidth', 2, 'DisplayName', 'G4 (4 alleles)');
    xlabel('Time', 'FontSize', 16);
    ylabel('Genotype Frequencies', 'FontSize', 16);
    title('Evolution of Genotype Frequencies', 'FontSize', 16);
    legend('FontSize', 14);
    set(gca, 'XScale', 'log', 'FontSize', 14);
    grid on;
    ylim([0 1]);
end

function dydt = odefunc(~, y, s, mu, nu, h1, h2, h3)
    % Extract log-transformed variables
    L00 = y(1);
    L01 = y(2);
    L11 = y(3);
    
    % Convert to original variables
    g00 = exp(L00);
    g01 = exp(L01);
    g11 = exp(L11);
    g10 = 1 - g00 - g01 - g11;  % constraint: g00 + g01 + g10 + g11 = 1
    
    % Mean fitness calculation
    wbar = 1 - 2*s*(h1*(g00*g10 + g00*g01) + h2*(g00*g11 + g01*g10) + h3*(g01*g11 + g10*g11)) - s*(h2*(g01^2 + g10^2) + g11^2);
    
    % Relative fitness calculations
    w0 = 1/wbar;
    w1 = (1 - s*h1)/wbar;
    w2 = (1 - s*h2)/wbar;
    w3 = (1 - s*h3)/wbar;
    w4 = (1 - s)/wbar;
    
    % Selection terms
    sel_g00 = g00^2*w0 + 0.5*2*g01*g00*w1 + 0.5*2*g10*g00*w1 + 0.25*(2*g00*g11 + 2*g01*g10)*w2;
    sel_g01 = 0.5*2*g01*g00*w1 + g01^2*w2 + 0.25*(2*g00*g11 + 2*g01*g10)*w2 + 0.5*2*g01*g11*w3;
    sel_g10 = 0.5*2*g10*g00*w1 + g10^2*w2 + 0.25*(2*g00*g11 + 2*g01*g10)*w2 + 0.5*2*g10*g11*w3;
    sel_g11 = 0.25*(2*g00*g11 + 2*g01*g10)*w2 + 0.5*2*g01*g11*w3 + 0.5*2*g10*g11*w3 + g11^2*w4;
    
    % Mutation terms
    mut_g00 = sel_g00*(1-mu)^2 + sel_g01*(1-mu)*nu + sel_g10*(1-mu)*nu + sel_g11*nu^2 - g00;
    mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu)*(1-nu) + sel_g10*mu*nu + sel_g11*(1-nu)*nu - g01;
    mut_g11 = sel_g00*mu^2 + sel_g01*mu*(1-nu) + sel_g10*mu*(1-nu) + sel_g11*(1-nu)^2 - g11;
    
    % mut_g10 is determined by constraint: mut_g10 = -mut_g00 - mut_g01 - mut_g11
    % This maintains the constraint that g00 + g01 + g10 + g11 = 1
    % mut_g10 = -mut_g00 - mut_g01 - mut_g11;
    
    % Log-transformed ODEs
    dL00dt = (1 / g00) * mut_g00;
    dL01dt = (1 / g01) * mut_g01;
    dL11dt = (1 / g11) * mut_g11;
    
    dydt = [dL00dt; dL01dt; dL11dt];
end

% Run the simulation
solve_allo_ode