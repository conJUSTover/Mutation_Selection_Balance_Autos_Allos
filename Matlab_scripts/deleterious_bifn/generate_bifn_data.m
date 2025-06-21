s_val_range = logspace(-9, -3, 1000);
mu_val = 2e-8;
nu_val = 1e-9;

% recessive case
h_val = 0;
h1_val = 0;
h2_val = 0;
h3_val = 0;

[dip_rec] = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
disp('Diploid data complete')

[auto_rec] = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Auto data complete')

[allo_rec] = allo_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Allo data complete')

writematrix(dip_rec, 'dip_rec.csv')
writematrix(auto_rec, 'auto_rec.csv')
writematrix(allo_rec, 'allo_rec.csv')

disp('Recessive case complete')

% additive case
h_val = .5;
h1_val = .25;
h2_val = .5;
h3_val = .75;

[dip_add] = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
disp('Diploid data complete')

[auto_add] = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Auto data complete')

[allo_add] = allo_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Allo data complete')

writematrix(dip_add, 'dip_add.csv')
writematrix(auto_add, 'auto_add.csv')
writematrix(allo_add, 'allo_add.csv')


disp('Additive case complete')

% dominant case
h_val = 1;
h1_val = 1;
h2_val = 1;
h3_val = 1;

[dip_dom_neutral, dip_dom_selected, dip_dom_unstable] = diploids_bifn_data(s_val_range, mu_val, nu_val, h_val);
disp('Diploid data complete')

[auto_dom_neutral, auto_dom_selected, auto_dom_unstable] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Auto data complete')

[allo_dom_neutral, allo_dom_selected, allo_dom_unstable] = allo_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp('Allo data complete')

writematrix(dip_dom_neutral, 'dip_dom_neutral.csv')
writematrix(dip_dom_selected, 'dip_dom_selected.csv')
writematrix(dip_dom_unstable, 'dip_dom_unstable.csv')

writematrix(auto_dom_neutral, 'auto_dom_neutral.csv')
writematrix(auto_dom_selected, 'auto_dom_selected.csv')
writematrix(auto_dom_unstable, 'auto_dom_unstable.csv')

writematrix(allo_dom_neutral, 'allo_dom_neutral.csv')
writematrix(allo_dom_selected, 'allo_dom_selected.csv')
writematrix(allo_dom_unstable, 'allo_dom_unstable.csv')

