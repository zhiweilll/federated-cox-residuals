test_25 = run_single_simulation_scenario(scenario_id = 25, seed = 12, n_total = 10000, n_sites = 3,
                               intercept_range = c(0, 0.5),
                               beta_fixed_range = c(0,0),
                               beta_tde_range = c(1, 1),
                               basehaz_range = c(0.01, 0.05),
                               tdefun_library = get_tdefun_library()) 


test_2 = run_single_simulation_scenario(scenario_id = 25, seed = 13, n_total = 10000, n_sites = 3,
                                         intercept_range = c(0.4, 0.6),
                                         beta_fixed_range = c(0,0),
                                         beta_tde_range = c(1.5, 2.5),
                                         basehaz_range = c(0.01, 0.05),
                                         tdefun_library = get_tdefun_library()) 


test_25 = replicate_simulation_scenario(
  scenario_id = 25,
  replicate_seed = 12,
  n_total = 10000,
  n_sites = 2,
  n_rep = 100,
  intercept_range = c(0, 0.2),
  beta_fixed_range = c(0,0),
  beta_tde_range = c(1, 1),
  basehaz_range = c(0.01, 0.05),
  tdefun_library = get_tdefun_library()) 


replicate_25 = replicate_simulation_scenario(
  scenario_id = 25,
  replicate_seed = 12,
  n_total = 10000,
  n_sites = 2,
  n_rep = 10,
  intercept_range = c(0, 0.2),
  beta_fixed_range = c(0,0),
  beta_tde_range = c(1, 1),
  basehaz_range = c(0.01, 0.05),
  tdefun_library = get_tdefun_library()) 



replicate_25_2 = replicate_simulation_scenario(
  scenario_id = 25,
  replicate_seed = 12,
  n_total = 10000,
  n_sites = 2,
  n_rep = 20,
  intercept_range = c(0, 0.2),
  beta_fixed_range = c(0,0),
  beta_tde_range = c(1, 1),
  basehaz_range = c(0.01, 0.05),
  tdefun_library = get_tdefun_library()) 


replicate_25_2$summary





test_25_2 = run_single_simulation_sites_first(scenario_id = 25, seed = 13, n_sites = 2,n_per_site = 4000,
                                  intercept_range = c(0.4, 0.6),
                                  beta_fixed_range = c(0,0),
                                  beta_tde_range = c(1.5, 2.5),
                                  basehaz_range = c(0.01, 0.05),
                                  tdefun_library = get_tdefun_library()) 
  