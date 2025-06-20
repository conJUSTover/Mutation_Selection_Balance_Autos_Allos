function stable_data = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val)
%Generates diploid bifurcation plotting data
%   for diploids, a nonlinear model and analysis
%   returns:
%       stable_data: matrix of s, q, g0, g1, and w-bar data for diploid
%                      equilibria which are stable


g0 = sym('g0');
g1 = sym('g1');
G0 = sym('G0');
G1 = sym('G1');
G2 = sym('G2');
s = sym('s');
h = sym('h');
mu = sym('mu');
nu = sym('nu');

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1 & 'Real');
assume(g1>=0 & g1<=1);
assume(s>=-1 & s<=1);
assume(h>=0 & h<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);

%%% g0 is p and g1 is q

% equations to parameterize relative fitnesses
w_bar = 1 - s*(h*2*g0*g1 + g1^2);
w0 = 1/w_bar;
w1 = (1-s*h)/w_bar;
w2 = (1-s)/w_bar;

% equations for selection
sel_g0 = w1*g1*g0 + w0*g0^2; 
sel_g1 = w2*g1^2 + w1*g1*g0;

% equations for mutation
mut_g0 = sel_g0*(1-mu) + sel_g1*nu - g0;
mut_g1 = sel_g0*mu + sel_g1*(1-nu) - g1;

%removing g1 from the equations
mut_g0 = subs(mut_g0, g1, 1-g0);

%derivative for linear stability analysis
g0_deriv = diff(mut_g0, g0);

[mut_g0_eval] = diff_eq_eval(mut_g0, mu, mu_val, nu, nu_val, h, h_val);

[g0_deriv_eval] = derivative_eval(g0_deriv, mu, mu_val, nu, nu_val, h, h_val);

stable_g0 = zeros(1, length(s_val_range));
stable_s = s_val_range;

for i = 1:length(s_val_range)

    [g0_roots] = dip_root_solns(mut_g0_eval, s, s_val_range(i), g0);

    [root_stabilities] = dip_linear_stability_analysis(g0_deriv_eval, s, s_val_range(i), g0, g0_roots);

    for j = 1:length(root_stabilities)
        if root_stabilities(j) == 1
            stable_g0(1, i) = g0_roots(j);
        end
    end
end

stable_g1 = ones(1, length(stable_g0)) - stable_g0;

stable_genotypes = [stable_g0.^2; 2*stable_g0.*stable_g1; stable_g1.^2];

stable_abs_fitness = [ones(1, length(stable_s)); ones(1, length(stable_s)) - h_val*stable_s; ones(1, length(stable_s)) - stable_s];

stable_w_bar = zeros(1, length(stable_s));

for i = 1:length(stable_w_bar)
    stable_w_bar(i) = dot(stable_genotypes(:, i), stable_abs_fitness(:, i));
end

stable_q = stable_g1;

stable_data_transposed = vertcat(stable_s, stable_q, stable_g0, stable_g1, stable_w_bar);
stable_data = stable_data_transposed.';

end

%%%%%%%%%%%
%FUNCTIONS%
%%%%%%%%%%%

function [g0_eval] = diff_eq_eval(mut_eqn_g0, mu, mu_val, nu, nu_val, h, h_val)

    g0_eval = subs(mut_eqn_g0, mu, mu_val);
    g0_eval = subs(g0_eval, nu, nu_val);
    g0_eval = subs(g0_eval, h, h_val);

end

function [g0_root] = dip_root_solns(g0_eqn_eval, s, s_val, g0)

    g0_eqn_eval = subs(g0_eqn_eval, s, s_val);
    
    g0_root = vpasolve(g0_eqn_eval, g0);
end

function [g0_derivative] = derivative_eval(g0_derivative, mu, mu_val, nu, nu_val, h, h_val)

    % creates a partially evaluated derivative for the constant parameters
    % (i.e. subs mu, nu, and h into the derivative since these are
    % constant parameters)

    g0_derivative = subs(g0_derivative, h, h_val);
    g0_derivative = subs(g0_derivative, mu, mu_val);
    g0_derivative = subs(g0_derivative, nu, nu_val);

end

function [fixed_pt_stabilities] = dip_linear_stability_analysis(g0_derivative, s, s_val, g0, g0_roots)

    g0_derivative = subs(g0_derivative, s, s_val);
    
    fixed_pt_stabilities = g0_roots;

    for i = 1:length(fixed_pt_stabilities)
        g0_derivative_eval = subs(g0_derivative, g0, g0_roots(i));
        if g0_derivative_eval < 0
            fixed_pt_stabilities(i) = 1; %1 indicates a stable node
        elseif g0_derivative_eval > 0
            fixed_pt_stabilities(i) = 0; %0 indicates an unstable node
        else
            disp('Error. Linear stabilitiy analysis failed. Derivative is 0.')
        end
    end


end