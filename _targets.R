library(tarchetypes)
library(tidyverse)
library(targets)

source("R/esm_helper.R")

library(crew)
tar_option_set(
  controller = crew_controller_local(workers = 8)
)

N_SIM = 100
SITES = c("Boon", "Craw", "Sali", "Merc", "Chri", "Ogle")

tar_option_set(packages = c("tidyverse"))

data_values <- bind_rows(
  expand_grid(type="train", site=SITES, depths=list(c(30, 35))),
  expand_grid(type="train", site=SITES, depths=list(c(30, 60))),
  expand_grid(type="train", site=SITES, depths=list(c(30))),
  expand_grid(type="train", site=SITES, depths=list(c(60))),
  expand_grid(type="test", site=SITES, depths=list(c(15,25,27.5,30,32.5,35,45))),
  expand_grid(type="test", site=SITES, depths=list(30+c(15,25,27.5,30,32.5,35,45))),
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

train_targets = list(
  # measure each section with noise
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
    tar_target(train_compAll_noise,
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
    mutate(train_name="train_%s%s_%s",
           test_name="test_%s_%s") %>%
    expand_grid(noise = "noise_",
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
  
  tibble(site=SITES) %>%
    mutate(train_name="train_comp_%s%s_%s",
           test_name="test_%s_%s") %>%
    expand_grid(noise = "noise_",
                train_depths = c("30_35"),
                test_depths = "15_25_27.5_30_32.5_35_45") %>%
    transmute(train_name = sprintf(train_name, noise, site, train_depths),
              test_name = sprintf(test_name, site, test_depths)),
  
  tibble(site=SITES) %>%
    mutate(train_name="train_compAll_%s%s_%s",
           test_name="test_%s_%s") %>%
    expand_grid(noise = "noise_",
                train_depths = c("30_60"),
                test_depths = c("15_25_27.5_30_32.5_35_45",
                                "45_55_57.5_60_62.5_65_75") ) %>%
    transmute(train_name = sprintf(train_name, noise, site, train_depths),
              test_name = sprintf(test_name, site, test_depths)),
  ) %>%
  mutate(train=rlang::syms(train_name),
         test=rlang::syms(test_name),
         truth=rlang::syms(str_replace(test_name, "test", "truth")))

predict_targets = tar_map(
  values=predict_values,
  names=c(train, test),
  
  # hyman
  tar_target(fit_hyman,
             train %>% 
               ensure_sim %>%
               nest_by(location_id, .sim) %>%
               mutate(fit=list(splinefun(data$Mineral_cumsum, data$SOCd_cumsum, method="hyman")))
  ),
  tar_target(predict_hyman,
             predict_fun(fit_hyman, test)),

  # linear but always use bottom layer
  tar_target(fit_linear2,
             train %>% 
               ensure_sim %>%
               nest_by(location_id, .sim) %>%
               mutate(fit=list(approxExtrapFun(tail(data$Mineral_cumsum,2), 
                                               tail(data$SOCd_cumsum, 2)
                                         )))),
  tar_target(predict_linear2,
             predict_fun(fit_linear2, test)),
  # fixed depth
  tar_target(fit_fixed,
            train %>%
                ensure_sim %>%
                filter(Mineral_cumsum > 0) %>% # fixed takes nearest value but never want zero
                nest_by(location_id, .sim) %>%
                mutate(fit=list(fixedfun(data$Mineral_cumsum, data$SOCd_cumsum)))),
  tar_target(predict_fixed,
             predict_fun(fit_fixed, test)),
  tar_target(predict,
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
    predict_targets$predict,
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
  )
)