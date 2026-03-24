library(targets)
library(tarchetypes)
library(tidyverse)

source("R/esm_helper.R")
source("R/targets_helper.R")
source("R/resample_helper.R")

library(crew)
tar_option_set(
  controller = crew_controller_local(workers = 4)
)

N_SIM = 100
RESAMPLE_SITES = c("Boon", "Craw", "Sali", "Merc", "Chri", "Ogle")

tar_option_set(packages = c("tidyverse", "sf"))

data_values <- bind_rows(
  expand_grid(type="train", site=RESAMPLE_SITES, depths=list(c(30, 35))),
  expand_grid(type="train", site=RESAMPLE_SITES, depths=list(c(30, 60))),
  expand_grid(type="train", site=RESAMPLE_SITES, depths=list(c(30))),
  expand_grid(type="test", site=RESAMPLE_SITES, depths=list(c(15,25,27.5,30,32.5,35,45))),
  
  expand_grid(type="train", site=RESAMPLE_SITES, depths=list(c(60))),
  expand_grid(type="test", site=RESAMPLE_SITES, depths=list(30+c(15,25,27.5,30,32.5,35,45))),
) %>% 
  rowwise() %>%
  mutate(depths_name = paste0(depths, collapse="_")) %>%
  mutate(measurements = c(rlang::sym(paste0("measurements_", site))))

truth_targets = tar_map(
  data_values %>% filter(type=="test"),
  names=c(site, depths_name),
  tar_target(test,
             measurements %>%
               filter(sample_depth_max %in% c(depths)) %>%
               select(location_id, Mineral_cumsum)),
  tar_target(truth,
             measurements %>%
               filter(sample_depth_max %in% c(depths)))
)

# "old" interpolation models (splines and linear)
# this sets up which train/test sets to run them on
old_data_values <- bind_rows(
  tibble(site=RESAMPLE_SITES) %>%
    mutate(train_name="train_%s%s_%s",
           test_name="test_%s_%s") %>%
    expand_grid(noise = c("", "noise_", "noise2_")[[2]],
                bind_rows(
                  expand_grid(
                    train_depths = c("30", "30_35", "30_60"),
                    test_depths = "15_25_27.5_30_32.5_35_45"),
                  expand_grid(
                    train_depths = c("60", "30_60"),
                    test_depths = "45_55_57.5_60_62.5_65_75")
                  )) %>%
    transmute(train_name = sprintf(train_name, noise, site, train_depths),
              test_name = sprintf(test_name, site, test_depths)),
  
  tibble(site=RESAMPLE_SITES) %>%
    mutate(train_name="train_comp_%s%s_%s",
           test_name="test_%s_%s") %>%
    expand_grid(noise = c("", "noise_", "noise2_")[[2]],
                train_depths = c("30_35"),
                test_depths = "15_25_27.5_30_32.5_35_45") %>%
    transmute(train_name = sprintf(train_name, noise, site, train_depths),
              test_name = sprintf(test_name, site, test_depths)),
  
  tibble(site=RESAMPLE_SITES) %>%
    mutate(train_name="train_compAll_%s%s_%s",
           test_name="test_%s_%s") %>%
    expand_grid(noise = c("", "noise_", "noise2_")[[2]],
                train_depths = c("30_60"),
                test_depths = c("15_25_27.5_30_32.5_35_45",
                                "45_55_57.5_60_62.5_65_75") ) %>%
    transmute(train_name = sprintf(train_name, noise, site, train_depths),
              test_name = sprintf(test_name, site, test_depths)),
  ) %>%
  mutate(train=rlang::syms(train_name),
         test=rlang::syms(test_name),
         truth=rlang::syms(str_replace(test_name, "test", "truth")))

old_targets = tar_map(
  values=old_data_values,
  names=c(train, test),
  tar_target(
    test_minus_train,
    test#get_test_minus_train(test, train)
  ),
  
  # hyman
  tar_target(fit_hyman,
             train %>% 
               ensure_sim %>%
               nest_by(location_id, .sim) %>%
               mutate(fit=list(splinefun(data$Mineral_cumsum, data$SOCd_cumsum, method="hyman")))
  ),
  tar_target(predict_hyman,
             predict_fun(fit_hyman, test_minus_train)),

  # linear but always use correction layer
  tar_target(fit_linear2,
             train %>% 
               ensure_sim %>%
               nest_by(location_id, .sim) %>%
               mutate(fit=list(approxExtrapFun(tail(data$Mineral_cumsum,2), 
                                               tail(data$SOCd_cumsum, 2)
                                         )))),
  tar_target(predict_linear2,
             predict_fun(fit_linear2, test_minus_train)),
  # fixed depth
  tar_target(fit_fixed,
            train %>%
                ensure_sim %>%
                filter(Mineral_cumsum > 0) %>% # fixed takes nearest value but never want zero
                nest_by(location_id, .sim) %>%
                mutate(fit=list(fixedfun(data$Mineral_cumsum, data$SOCd_cumsum)))),
  tar_target(predict_fixed,
             predict_fun(fit_fixed, test_minus_train)),
  tar_target(predict_old,
             bind_rows(`hyman`=predict_hyman,
                       `linear2`=predict_linear2,
                       `fixed`=predict_fixed,
                       .id="name") %>%
               mutate(train=train_name, test=test_name) %>%
               group_by(name, location_id) %>%
               arrange(Mineral_cumsum) %>%
               # add marginal SOC prediction
               compute_marginal_pred
  )
)



list(
  tar_target(resample_soils,
             get_resample_soils()),
  tar_target(resample_labs,
             get_resample_labs()),
  tar_target(measurements_resample,
             get_resample_measurements(resample_labs, resample_soils)
             ),
  tar_map(values=tibble(site=RESAMPLE_SITES),
          tar_target(measurements,
                     measurements_resample %>%
                       filter(str_starts(location_id, site)) %>%
                       #filter(sample_depth_max <= 45) %>%
                       group_by(location_id) %>%
                       #filter(n() == 7) %>%
                       ungroup %>%
                       add_ESM_columns)
          ),
  
  tar_map(
    data_values %>% filter(type=="train"),
    names=c(site, depths_name),
    tar_target(train,
               get_train(measurements, depths)),
    tar_target(train_noise,
               lapply(1:N_SIM, function(i) 
                 train %>% 
                   add_noise2(0.04, 0.02) %>% 
                   mutate(.sim = i)) %>%
                 do.call(bind_rows, .))
  ),
  
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
  
  tar_map(
    data_values %>% filter(type=="train", depths_name %in% c("30_60")),
    names=c(site, depths_name),
    tar_target(train_compAll_noise,
               lapply(1:N_SIM, function(i) get_train(measurements, depths) %>%
                        add_noise2(0.04, 0.02) %>%
                        pool_all %>%
                        mutate(.sim = i)) %>%
                 do.call(bind_rows, .)
    )
  ),
    
  truth_targets,
  
  old_targets,
  tar_combine(
    predict_old_combined,
    old_targets$predict_old,
    command=bind_rows(!!!.x)
  ),
  
  tar_target(validation_combined,
             get_validation_combined(predict_old_combined,
                                     truth_combined)
  ),
  
  tar_combine(
    truth_combined,
    c(truth_targets$truth),
    command=bind_rows_with_name(!!!.x)
  )
)
