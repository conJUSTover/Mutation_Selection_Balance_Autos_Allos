% for autos, a visualization of the phase plane for different s values
% shows basins of attraction and how bifurcations play out

%%% Recreates Figure S5

s_val_range = [1e-9, 3e-7, .3]; % 3 element vector of s_vals

mu_val = 2e-8; % constant value of forward mutation rate
nu_val = 1e-9; % constant value of backward mutation rate
mut_ratio_val = mu_val/nu_val; % ratio of forward to backward mutation rate
a_val = 0; % constant value of alpha (double reduction rate)

h1_val = 1; % h1 dominance coefficient value, constant
h2_val = 1; % h2 dominance coefficient value, constant
h3_val = 1; % h3 dominance coefficient value, constant

syms a s q G0 G1 G2 G3 G4 g0 g1 g2 h1 h2 h3 mu nu

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(g2>=0 & g2<=1);
assume(s>=-1 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(a>=0 & a<=1/6);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);
assume(G3>=0 & G3<=1);
assume(G4>=0 & G4<=1);

% equations to parameterize relative fitnesses
wbar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4);
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g0 = G0*w0+(1/2 + a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (a/4)*G3*w3;
sel_g1 = (1/2 - a/2)*G1*w1 + (2/3 - 2*a/3)*G2*w2 + (1/2 - a/2)*G3*w3;
sel_g2 = (a/4)*G1*w1 + (1/6 + a/3)*G2*w2 + (1/2 + a/4)*G3*w3 + G4*w4;

% equations for mutation
mut_g0 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2) - g0;
mut_g1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu - g1;
mut_g2 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2) - g2;

mut_exp_set = [mut_g0, mut_g1, mut_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(mut_exp_set)
    mut_exp_set(i) = subs(mut_exp_set(i), G0, g0^2);
    mut_exp_set(i) = subs(mut_exp_set(i), G1, 2*g0*g1);
    mut_exp_set(i) = subs(mut_exp_set(i), G2, (2*g0*g2 + g1^2));
    mut_exp_set(i) = subs(mut_exp_set(i), G3, 2*g1*g2);
    mut_exp_set(i) = subs(mut_exp_set(i), G4, g2^2);
    mut_exp_set(i) = subs(mut_exp_set(i), g2, (1-g1-g0));
end

%creates the Jacobian of the system
jacobian_1 = [diff(mut_exp_set(1), g0), diff(mut_exp_set(1), g1); diff(mut_exp_set(2), g0), diff(mut_exp_set(2), g1)];


figure;

tiledlayout(1, 3, "TileSpacing", "compact");

for h = 1:3
    % Move to the next tile for each iteration
    nexttile;
    
    if h == 1
        disp('s<mu')
    elseif h == 2
        disp('s>mu')
    else
        disp('s>>mu')
    end

    %solves for the fixed points of the system
    [g0_value, g1_value] = numeric_solver(mut_exp_set(1), mut_exp_set(2), mu, mu_val, nu, nu_val, s, s_val_range(h), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g1);

    %initializes equations to create vectors for the quiver plot
    [g0_eqn, g1_eqn] = quiver_plot_init(mut_exp_set(1), mut_exp_set(2), mu, mu_val, nu, nu_val, s, s_val_range(h), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val);

    x_1 = linspace(0, 1, 100);
    y_1 = zeros(1, 100);
    y_2 = zeros(1, 100);

    for i = 1:length(x_1)
         y_1_soln = vpasolve(subs(g0_eqn, g0, x_1(i)), g1);
         for j = 1:length(y_1_soln)    
             if y_1_soln(j) >= 0 && y_1_soln(j) <= 1 && imag(y_1_soln(j)) == 0 
                 y_1(i) = y_1_soln(j);
             end
         end
         y_2_soln = vpasolve(subs(g1_eqn, g0, x_1(i)), g1);
         for j = 1:length(y_2_soln)    
             if y_2_soln(j) >= 0 && y_2_soln(j) <= 1 && imag(y_2_soln(j)) == 0 
                 y_2(i) = y_2_soln(j);
             end
         end
    end

    null_1 = plot(x_1, y_1, 'LineWidth', 7.5);
    null_1.Color = "#2a81a5";
    hold on
    null_2 = plot(x_1, y_2, 'LineWidth', 4);
    null_2.Color = "#e57c34";

    y_3 = 1-x_1;

    plot(x_1, y_3, 'LineWidth', 2, 'Color', [.5, .5, .5])

    y_HWE = 2*(x_1.^.5 - x_1);
    plot(x_1, y_HWE, 'LineWidth', 1, 'Color', 'k')

    for i = 1:length(g0_value)

        jacobian_eval = zeros(length(jacobian_1));
        %evaluating the jacobian
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
            jacobian_eval(j, k) = pd_evaluation(jacobian_1(j, k), mu, mu_val, nu, nu_val, s, s_val_range(h), h1, h1_val, h2, h2_val, h3, h3_val, a, a_val, g0, g0_value(i), g1, g1_value(i)); 
            end
        end

        det_jac = det(jacobian_eval);
        if det_jac < 0
            plot(g0_value(i), g1_value(i), 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerSize', 12)
            [eigenvectors, eigenvalues] = eig(jacobian_eval);
            eigenvalues = sum(eigenvalues);
            for k = 1:2
                if eigenvalues(k) < 0
                    x = linspace(g0_value(i)-.4, g0_value(i)+.4, 10);
                    slope = eigenvectors(2, k)/eigenvectors(1, k);
                    y = slope*(x - g0_value(i)) + g1_value(i);
                    plot(x, y, "LineStyle","--", "Color", 'k', 'LineWidth', 2)
                end
            end
        else
            plot(g0_value(i), g1_value(i), 'o', 'Color', 'k', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 12)
        end     

        %creates a grid of values covering [0, 1] x [0, 1]
        g0_input_values = 0:1/(iterations-1):1;
        g1_input_values = 0:1/(iterations-1):1;

        %creates coordinate data for the quiver plot using meshgrid
        [g0_coordinates, g1_coordinates] = meshgrid(g0_input_values, g1_input_values);

        %creates a blank array to store vector values for the quiver plot
        g0_vector_values = zeros(iterations, iterations);
        g1_vector_values = zeros(iterations, iterations);

        %generates the vectors for the quiver plot
        for j = 1:iterations
            for k = 1:iterations
                [g0_vector, g1_vector] = quiver_plot_vectors(g0_eqn, g1_eqn, g0_coordinates(j, k), g1_coordinates(j, k), g0, g1);
                %restricting the output to contain only the space where the
                %sum of gamete frequencies is less than 1
                if g0_coordinates(j, k) + g1_coordinates(j, k) <= 1
                    g0_vector_values(j, k) = g0_vector;
                    g1_vector_values(j, k) = g1_vector;
                else
                    g0_vector_values(j, k) = 0;
                    g1_vector_values(j, k) = 0;
                end
            end
        end

        if isscalar(g0_value)
            quiver(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, 'Color', [.5 .5 .5])
        elseif i == 3
            quiver(g0_coordinates, g1_coordinates, g0_vector_values, g1_vector_values, 'Color', [.5 .5 .5])
        end

    end

    xlim([0 1])
    ylim([0 1])

    if h == 1
        ylabel('g_1')
        yticks([0, .25, .5, .75, 1])
        xticks([0, .25, .5, .75, 1])
        text(gca, .98, .98, 'A', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 20);
    end 

    if h == 2
        lgd = legend({'g_0 Nullcline', 'g_1 Nullcline', '', 'HWE', 'Stable Node', '', 'Saddle Point', 'Eigenvector', 'Vector Field'});
        xlabel('g_0')
        set(gca, 'YTickLabel', [])
        xticks([0, .25, .5, .75, 1])
        text(gca, .98, .98, 'B', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 20);
    end

    if h == 3
        set(gca, 'YTickLabel', [])
        xticks([0, .25, .5, .75, 1])
        text(gca, .98, .98, 'C', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 20);
    end

    set(gca, 'FontSize', 18);
    set(gca, 'XTickLabelRotation', 0);
end
lgd.Layout.Tile = 'east';
sgtitle('Autotetraploid Phase Space Diagrams', 'FontSize', 24);
set(gcf, 'Units', 'inches');
set(gcf, 'Position', [0, 0, 13, 5]); % [left, bottom, width, height]
exportgraphics(gcf, 'auto_phase_space.pdf', 'ContentType','vector', 'Resolution',300)

% save figure to pdf
% run this from the command window to save the PDF after editing the figure
% in the interactive window
% exportgraphics(gcf, 'auto_phase_space.pdf', 'ContentType','vector', 'Resolution',300)


%function which uses vpasolve to evaluate the fixed points of the system
function [g0_value, g1_value] = numeric_solver(mut_g0_eqn, mut_g1_eqn, mu, mu_value, nu, nu_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g1)

    g0_eqn = subs(mut_g0_eqn, mu, mu_value);
    g0_eqn = subs(g0_eqn, nu, nu_value);
    g0_eqn = subs(g0_eqn, s, sel_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g1_eqn = subs(mut_g1_eqn, mu, mu_value);
    g1_eqn = subs(g1_eqn, nu, nu_value);
    g1_eqn = subs(g1_eqn, s, sel_value);
    g1_eqn = subs(g1_eqn, h1, h1_value);
    g1_eqn = subs(g1_eqn, h2, h2_value);
    g1_eqn = subs(g1_eqn, h3, h3_value);
    g1_eqn = subs(g1_eqn, a, a_value);


    [g0_value, g1_value] = vpasolve([g0_eqn, g1_eqn], [g0, g1]);

end

%function which acts as a partial derivative evaluation tool 
%used to evaluate the jacobian at each entry
function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_value, nu, nu_value, s, s_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g0_sub_value, g1, g1_sub_value)

    pd_value = subs(jacobian_entry, mu, mu_value);
    pd_value = subs(pd_value, nu, nu_value);
    pd_value = subs(pd_value, s, s_value);
    pd_value = subs(pd_value, h1, h1_value);
    pd_value = subs(pd_value, h2, h2_value);
    pd_value = subs(pd_value, h3, h3_value);
    pd_value = subs(pd_value, a, a_value);
    pd_value = subs(pd_value, g0, g0_sub_value);
    pd_value = subs(pd_value, g1, g1_sub_value);

end

%initializes the quiver plot by substituting in all of the parameter values
%which are constant (i.e. s, mu, h1, h2, h3, a)
function [g0_eqn, g1_eqn] = quiver_plot_init(mut_g0_eqn, mut_g1_eqn, mu, mu_value, nu, nu_value, s, sel_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value)

    g0_eqn = subs(mut_g0_eqn, mu, mu_value);
    g0_eqn = subs(g0_eqn, nu, nu_value);
    g0_eqn = subs(g0_eqn, s, sel_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g1_eqn = subs(mut_g1_eqn, mu, mu_value);
    g1_eqn = subs(g1_eqn, nu, nu_value);
    g1_eqn = subs(g1_eqn, s, sel_value);
    g1_eqn = subs(g1_eqn, h1, h1_value);
    g1_eqn = subs(g1_eqn, h2, h2_value);
    g1_eqn = subs(g1_eqn, h3, h3_value);
    g1_eqn = subs(g1_eqn, a, a_value);

end

%generates vectors for the quiver plot by substituting the current values
%of g0 and g1
function [g0_vector, g1_vector] = quiver_plot_vectors(g0_eqn, g1_eqn, g0_value, g1_value, g0, g1)

    g0_vector = subs(g0_eqn, g0, g0_value);
    g0_vector = subs(g0_vector, g1, g1_value);

    g1_vector = subs(g1_eqn, g0, g0_value);
    g1_vector = subs(g1_vector, g1, g1_value);

end

%for bifurcation analysis

function [pd_value] = bifn_pd_evaluation(jacobian_entry, mu, mu_value, nu, nu_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value)

    pd_value = subs(jacobian_entry, mu, mu_value);
    pd_value = subs(pd_value, nu, nu_value);
    pd_value = subs(pd_value, h1, h1_value);
    pd_value = subs(pd_value, h2, h2_value);
    pd_value = subs(pd_value, h3, h3_value);
    pd_value = subs(pd_value, a, a_value);

end

function [g0_bifn_value, g1_bifn_value, s_bifn_value] = bifn_numeric_solver(mut_g0_eqn, mut_g1_eqn, jacobian, mu, mu_value, nu, nu_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value, g0, g1, s, initial_conditions)

    g0_eqn = subs(mut_g0_eqn, mu, mu_value);
    g0_eqn = subs(g0_eqn, nu, nu_value);
    g0_eqn = subs(g0_eqn, h1, h1_value);
    g0_eqn = subs(g0_eqn, h2, h2_value);
    g0_eqn = subs(g0_eqn, h3, h3_value);
    g0_eqn = subs(g0_eqn, a, a_value);

    g1_eqn = subs(mut_g1_eqn, mu, mu_value);
    g1_eqn = subs(g1_eqn, nu, nu_value);
    g1_eqn = subs(g1_eqn, h1, h1_value);
    g1_eqn = subs(g1_eqn, h2, h2_value);
    g1_eqn = subs(g1_eqn, h3, h3_value);
    g1_eqn = subs(g1_eqn, a, a_value);

    jacobian_eval_bifn = jacobian;

    for j = 1:length(jacobian_eval_bifn)
        for k = 1:length(jacobian_eval_bifn)
            jacobian_eval_bifn(j, k) = bifn_pd_evaluation(jacobian(j, k), mu, mu_value, nu, nu_value, h1, h1_value, h2, h2_value, h3, h3_value, a, a_value); 
        end
    end

    [bifn_eigenvectors, bifn_eigenvalues] = eig(jacobian_eval_bifn);

    bifn_eig_det = det(bifn_eigenvalues) == 0;
    
    [g0_bifn_value, g1_bifn_value, s_bifn_value] = vpasolve([g0_eqn, g1_eqn, bifn_eig_det], [g0, g1, s], initial_conditions);

end

