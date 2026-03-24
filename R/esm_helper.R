# # make sure that there are no missing sections
# filter_complete_profiles = function(df) {
#   df %>%
#     group_by(location_id) %>%
#     mutate(adjacent=sample_depth_min == lag(sample_depth_max, default=0)) %>%
#     filter(cumprod(adjacent) == 1)
# }
# 
# Add SOCd, SOCd_cumsum, Mineral, Mineral_cumsum, and (0,0) row
add_ESM_columns = function(df) {
  df %>%
    mutate(SOCd = SOCc*BD*(sample_depth_max - sample_depth_min),
           Mineral = BD*(sample_depth_max - sample_depth_min)*100 - SOCd*1.724) %>%
    arrange(location_id, sample_depth_max) %>%
    group_by(location_id) %>%
    mutate(SOCd_cumsum = cumsum(SOCd),
           Mineral_cumsum = cumsum(Mineral)) %>%
    ungroup %>%
    bind_rows(tibble(location_id=unique(.$location_id),
                     sample_depth_max=0,
                     Mineral_cumsum=0,
                     SOCd_cumsum=0)) %>%
    arrange(location_id, sample_depth_max)
}
# 
# # if you filter out some rows, then the SOCd and Mineral rows are wrong, so this fixes them
# update_ESM_columns = function(df) {
#   df %>%
#     group_by(location_id) %>%
#     arrange(location_id, sample_depth_max) %>%
#     mutate(SOCd = SOCd_cumsum - lag(SOCd_cumsum, default=0),
#            Mineral = Mineral_cumsum - lag(Mineral_cumsum, default=0)) %>%
#     ungroup
# }
# 
# for making predictions using old ESM adjustment methods
predict_fun = function(nested_fits, test) {
  nested_fits %>%
    inner_join(test) %>%
    mutate(SOCd_cumsum_predict = fit(Mineral_cumsum)) %>%
    select(Mineral_cumsum, SOCd_cumsum_predict) %>%
    unnest()
}

# given a grouped, sorted df with Mineral_cumsum, SOCd_cumsum_predict
# calculate SOCd_predict
compute_marginal_pred = function(grouped_df) {
  grouped_df %>%
    mutate(SOCd_predict = (SOCd_cumsum_predict - lag(SOCd_cumsum_predict, default=0)))
}
# 
# compute_marginal_truth = function(grouped_df) {
#   grouped_df %>%
#     mutate(SOCd = (SOCd_cumsum - lag(SOCd_cumsum, default=0)))
# }
# 
closest_index <- function(x_new, x) {
  which.min(abs(x_new - x))
}

fixedfun = function(x, y) {
  env <- rlang::new_environment()
  env$x <- x
  env$y <- y

  myfunc = function(x_new_vec) {
    closest_indices = sapply(x_new_vec, closest_index, x = env$x)
    env$y[closest_indices]
  }

  myfunc
}

get_train = function(measurements, depths) {
  if(!"n_obs" %in% measurements) {
    measurements$n_obs = 1
  }

  measurements %>%
    left_join(tibble(sample_depth_max=c(0,depths),
                     depth=c(0,depths))) %>%
    group_by(location_id) %>%
    fill(depth, .direction="up") %>%
    filter(!is.na(depth)) %>% # drop deeper sections
    group_by(location_id, depth) %>%
    arrange(location_id, depth) %>%
    summarize(sample_depth_min = min(sample_depth_min),
              sample_depth_max = max(sample_depth_max),
              SOCd=sum(SOCd),
              SOCd_cumsum=last(SOCd_cumsum),
              n_obs=sum(n_obs),
              Mineral_cumsums = list(Mineral_cumsum),
              Mineral_cumsum=last(Mineral_cumsum)) %>%
    select(-depth) %>%
    # make sure all depths are present in the profile
    group_by(location_id) %>%
    filter(n() == length(c(0,depths))) %>%
    ungroup
}

add_noise2 = function(measurements, tau_SOC, tau_BD) {
  if(!"n_obs" %in% colnames(measurements)) {
    measurements$n_obs = 1
  }

  measurements %>%
    group_by(location_id) %>%
    mutate(Mineral = Mineral_cumsum - lag(Mineral_cumsum)) %>%
    ungroup %>%
    filter(sample_depth_max > 0) %>%
    mutate(SOCc = SOCd/(SOCd*1.724 + Mineral)*100,
           BD = SOCd/SOCc/(sample_depth_max - sample_depth_min)) %>%
    select(location_id, sample_depth_min, sample_depth_max, SOCc, BD, n_obs) %>%
    mutate(SOCc = SOCc *rnorm(n(), mean=1, sd=tau_SOC),
           BD = BD*rnorm(n(), mean=1, sd=tau_BD)) %>%
    add_ESM_columns()
}

rep_first_bottom = function(df) {
  df %>%
    mutate(bottom = sample_depth_max == max(sample_depth_max)) %>%
    group_by(sample_depth_max) %>%
    mutate(SOCc = ifelse(bottom, first(SOCc), SOCc),
           BD = ifelse(bottom, first(BD), BD)) %>%
    filter(sample_depth_max > 0 ) %>%
    add_ESM_columns()
}

# linear extrapolation in functional form (a combination of approxExtrap and approxfun)
approxExtrapFun = function(x,y) {
  env <- rlang::new_environment()
  env$x <- x
  env$y <- y

  function(newx) {
    t = (newx - env$x[[1]])/(env$x[[2]] - env$x[[1]])
    t*env$y[[2]] + (1-t)*env$y[[1]]
  }
}
# 
# Given ESM columns of SOCd and Mineral_cumsum
# add back the SOCc and BD
# requires sample_depth_max and sample_depth_min as well
add_SOCc_BD_to_ESM = function(ESM) {
  ESM %>%
    arrange(pick(any_of(".sim"), location_id, sample_depth_max)) %>%
    group_by(pick(any_of(".sim"), location_id)) %>%
    mutate(Mineral = Mineral_cumsum - lag(Mineral_cumsum)) %>%
    mutate(SOCc = SOCd / (Mineral + SOCd*1.724)*100,
           BD = (Mineral + SOCd*1.724)/(sample_depth_max - sample_depth_min)/100) %>%
    select(-Mineral) %>%
    ungroup
}
# 
# pool all measurements like secondary data analysis
# first average BD and SOCc
# then calculate mineral cumsum
# then replicate across locations
pool_all = function(ESM) {
  ESM %>%
    add_SOCc_BD_to_ESM() %>%
    group_by(sample_depth_min, sample_depth_max) %>%
    transmute(location_id,
              SOCc = mean(SOCc),
              BD = mean(BD)) %>%
    filter(sample_depth_min >= 0) %>% # deal with those initial rows...
    add_ESM_columns()
}

# pool the SOCc in the bottom layer
# note this code is weird, BD is not actually BD here... and SOCc is proportion not percent
# but it works so leaving it for now
pool_bottom = function(train) {
  train %>%
    group_by(location_id) %>%
    mutate(Mineral = Mineral_cumsum - lag(Mineral_cumsum)) %>%
    mutate(SOCc = SOCd / (Mineral + SOCd*1.724),
           BD = SOCd/SOCc) %>%
    group_by(sample_depth_max) %>%
    # Only observe mean SOC% and mea
    mutate(SOCc = weighted.mean(SOCc, BD),
           BD=mean(BD)) %>%
    ungroup %>%
    mutate(SOCd_pooled = SOCc*BD,#Mineral/(1-SOCc*1.724),
           .after=SOCd,
           SOCd = ifelse(sample_depth_max==max(sample_depth_max), SOCd_pooled, SOCd)) %>%
    select(-SOCc, -BD, SOCd_pooled) %>%
    group_by(location_id) %>%
    mutate(SOCd_cumsum=cumsum(coalesce(SOCd, 0))) %>%
    ungroup
}
# 
# if no .sim column make it equal to 1
ensure_sim = function(df) {
  if(!".sim" %in% colnames(df)) {
    df$.sim = 1
  }

  df
}
# 
get_validation_combined = function(predict_old_combined, truth_combined) {
  predict_old_combined %>%
    #bind_rows(tar_read(predict_gp_summary_combined)) %>%
    inner_join(truth_combined %>%
                 mutate(test=str_replace(name, "truth", "test")) %>%
                 select(-name)) %>%
    ungroup %>%
    filter(!grepl("minimax", train)) %>%
    extract(train, into=c("site", "train_depths"), regex="([A-Za-z]+)_([_0-9]+)$", remove =FALSE) %>%
    mutate(test_depths = str_extract(test, "[0-9_\\.]+$") %>% str_sub(2)) %>%
    mutate(n_depths = 1 + str_count(train_depths, "_"), .after="train_depths") %>%
    filter(site != "Alma") %>%
    # ignore fixed depth with correction layer
    filter(name != "fixed" | n_depths == 1) %>%
    # hyman only makes sense with correction layer
    filter(name != "hyman" | n_depths > 1) %>%
    filter(name != "linear") %>%
    #filter(name != "fixed" | sample_depth_max == 30) %>%
    mutate(error_cumsum=SOCd_cumsum_predict - SOCd_cumsum,
           error = SOCd_predict - SOCd
    )
}