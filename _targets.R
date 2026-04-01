library(tarchetypes)
library(tidyverse)
library(targets)

source("R/esm_helper.R")

library(crew)
tar_option_set(
  controller = crew_controller_local(workers = 8)
)

N_SIM = 200
SITES = c("Boon", "Craw", "Sali", "Merc", "Chri", "Ogle")

tar_option_set(packages = c("tidyverse"))

data_values <- bind_rows(
  tibble_row(type="train", depths=list(c(30, 35))),
  tibble_row(type="train", depths=list(c(30, 60))),
  tibble_row(type="train", depths=list(c(30))),
  tibble_row(type="train", depths=list(c(60))),
  tibble_row(type="test", depths=list(c(15,25,27.5,30,32.5,35,45))),
  tibble_row(type="test", depths=list(30+c(15,25,27.5,30,32.5,35,45))),
) %>% 
  expand_grid(site=SITES) %>%
  rowwise() %>%
  mutate(depths_name = paste0(depths, collapse="_")) %>%
  mutate(measurements = c(sym(paste0("measurements_", site))))

truth_targets = tar_map(
  data_values %>% filter(type=="test"),
  names=c(site, depths_name),
  # these are the test sets with location ID and Mineral cumsum
  tar_target(test,
             measurements %>%
               filter(sample_depth_max %in% c(depths)) %>%
               select(location_id, Mineral_cumsum)),
  # these are the ground truth sets including SOC stock
  tar_target(truth,
             measurements %>%
               filter(sample_depth_max %in% c(depths)))
)

train_targets = list(
  tar_map(
    data_values %>% filter(type=="train"),
    names=c(site, depths_name),
    # aggregate without measurement error
    tar_target(train,
               get_train(measurements, depths)),
    # add measurement error
    tar_target(train_ind_noise,
               lapply(1:N_SIM, function(i) 
                 train %>% 
                   add_noise2(0.04, 0.02) %>% 
                   mutate(.sim = i)) %>%
                 do.call(bind_rows, .))
  ),
  
  # composite correction layer (30-35 cm) in 0-30, 30-35 cm design
  tar_map(
    data_values %>% filter(type=="train", depths_name %in% c("30_35")),
    names=c(site, depths_name),
    tar_target(train_comp,
               get_train(measurements, depths) %>% pool_bottom),
    tar_target(train_comp_noise,
               lapply(1:N_SIM, function(i) train_comp %>%
                        add_noise2(0.04, 0.02) %>%
                        rep_first_bottom %>% 
                        mutate(.sim = i)) %>%
                 do.call(bind_rows, .)
    )
  ),
  
  # average both sections in 0-30, 30-60 cm design
  tar_map(
    data_values %>% filter(type=="train", depths_name %in% c("30_60")),
    names=c(site, depths_name),
    tar_target(train_avg_noise,
               lapply(1:N_SIM, function(i) get_train(measurements, depths) %>%
                        add_noise2(0.04, 0.02) %>%
                        pool_all %>%
                        mutate(.sim = i)) %>%
                 do.call(bind_rows, .)
    )
  )
)

# ESM methods (splines and linear)
# this sets up which train/test sets to run them on
predict_values <- bind_rows(
  tibble(site=SITES) %>%
    expand_grid(bind_rows(
                  expand_grid(
                    train_depths = c("30", "30_35", "30_60"),
                    test_depths = "15_25_27.5_30_32.5_35_45"),
                  expand_grid(
                    train_depths = c("60", "30_60"),
                    test_depths = "45_55_57.5_60_62.5_65_75")
                  )) %>%
    mutate(measurement_type="ind"),
  
  tibble(site=SITES) %>%
    expand_grid(train_depths = c("30_35"),
                test_depths = "15_25_27.5_30_32.5_35_45") %>%
    mutate(measurement_type="comp"),
  
  tibble(site=SITES) %>%
    expand_grid(train_depths = c("30_60"),
                test_depths = c("15_25_27.5_30_32.5_35_45",
                                "45_55_57.5_60_62.5_65_75") ) %>%
    mutate(measurement_type="avg"),
  ) %>%
  transmute(
    measurement_type,
    train_name = str_glue("train_{measurement_type}_noise_{site}_{train_depths}"),
    test_name = str_glue("test_{site}_{test_depths}"),
    train_n_depths = str_count(train_depths, "_") + 1) %>%
  mutate(train=syms(train_name),
         test=syms(test_name),
         truth=syms(str_replace(test_name, "test", "truth")))

predict_values = predict_values %>% 
   expand_grid(tibble(
     method = c("linear2", "hyman", "fixed"),
     min_train_depths = c(1, 2, 1),
     max_train_depths = c(Inf, Inf, 1))) %>%
   filter(train_n_depths >= min_train_depths,
          train_n_depths <= max_train_depths) %>%
  mutate(fit_function = syms(str_glue("fit_{method}")))

predict_targets = tar_map(
  values=predict_values,
  names=c(method, train, test),
  
  tar_target(fit_long,
             train %>% 
               nest_by(location_id, .sim) %>%
               mutate(fit=list(fit_function(data)))
  ),
  tar_target(predict_long,
             predict_fun(fit_long, test)),
  
  # add provenance and marginal SOC predictions
  tar_target(predict_extra,
             predict_long%>%
               mutate(name=method, 
                      train=train_name, 
                      test=test_name,
                      measurement_type=measurement_type) %>%
               group_by(location_id) %>%
               arrange(Mineral_cumsum) %>%
               # add marginal SOC prediction
               compute_marginal_pred)
)

list(
  tar_target(measurements_file,
             "data/measurements.csv", format="file"),
  tar_target(measurements_all,
             read_csv(measurements_file)
             ),
  tar_map(values=tibble(site=SITES),
          tar_target(measurements,
                     measurements_all %>%
                       filter(str_starts(location_id, site)) %>%
                       group_by(location_id) %>%
                       ungroup %>%
                       add_ESM_columns)
          ),
    
  train_targets,
  truth_targets,
  
  predict_targets,
  
  tar_combine(
    predict_combined,
    predict_targets$predict_extra,
    command=bind_rows(!!!.x)
  ),
  
  tar_combine(
    truth_combined,
    c(truth_targets$truth),
    command=bind_rows_with_name(!!!.x)
  ),
  
  tar_target(validation_combined,
             get_validation_combined(predict_combined,
                                     truth_combined)
  ),
  tar_target(scenarios,
             get_scenarios()),
  tar_target(
    validation_combined_scenario,
    scenarios %>%
      expand_grid(.sim = 1:N_SIM, core_id=1:12) %>%
      group_by(scenario, target_depth, .sim, core_id) %>%
      sample_n(1, weight=weight) %>%
      select(-weight) %>%
      inner_join(validation_combined %>%
                   mutate(target_depth = ifelse(grepl("30", test), 30, 60)) %>%
                   mutate(core_id = as.numeric(str_split_i(location_id, "_", 2))))
  )
)