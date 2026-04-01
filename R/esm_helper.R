bind_rows_with_name = function(..., .names_to="name") {
  args = list(...)
  dots = match.call(expand.dots = FALSE)$...
  names = sapply(dots, deparse)
  for(i in 1:length(args)) {
    args[[i]][.names_to] = names[[i]]
  }
  df = bind_rows(args)
  df
}

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

# for making predictions using old ESM adjustment methods
predict_fun = function(nested_fits, test) {
  nested_fits %>%
    inner_join(test, relationship="many-to-many") %>%
    mutate(SOCd_cumsum_predict = fit(Mineral_cumsum)) %>%
    select(Mineral_cumsum, SOCd_cumsum_predict) %>%
    unnest(c())
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

fit_hyman = function(data) {
  splinefun(data$Mineral_cumsum, data$SOCd_cumsum, method="hyman")
}

fit_fixed = function(data) {
  data = data %>% filter(Mineral_cumsum > 0)
  fixedfun(data$Mineral_cumsum, data$SOCd_cumsum)
}

fit_linear2 = function(data) {
  approxExtrapFun(tail(data$Mineral_cumsum,2), 
                  tail(data$SOCd_cumsum, 2))
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

get_validation_combined = function(predict_combined, truth_combined) {
  predict_combined %>%
    inner_join(truth_combined %>%
                 mutate(test=str_replace(name, "truth", "test")) %>%
                 select(-name)) %>%
    ungroup %>%
    extract(train, into=c("site", "train_depths"), regex="([A-Za-z]+)_([_0-9]+)$", remove =FALSE) %>%
    mutate(test_depths = str_extract(test, "[0-9_\\.]+$") %>% str_sub(2)) %>%
    mutate(error_cumsum=SOCd_cumsum_predict - SOCd_cumsum,
           error = SOCd_predict - SOCd
    )
}

discretize_bd_dist = function(slices, mean, sd_log) {
  targets = c(slices$sample_depth_max)
  mids = (targets[-1] + targets[-length(targets)]) / 2
  bounds = c(0, mids, Inf)
  probs = diff(plnorm(bounds, mean = log(mean) - sd_log^2/2, sd = sd_log))
  tibble(depth=targets, weight=probs)
}
 
get_scenarios = function() {
  slices = tibble(sample_depth_max = c(15,25,27.5,30,32.5,35,45,55,57.5,60,62.5,65,75)) %>%
    mutate(sample_depth_min = lag(sample_depth_max, default=0)) %>%
    mutate(xmin = (sample_depth_max + sample_depth_min)/2,
           xmax = coalesce(lead(xmin), 80))
  
 scenarios0 = expand_grid(tibble(
    scenario=c("Compact", "Expand", "More Compact", "More Expand"), 
    dd = c(-2.5, +2.5, -5, +5) ),
    tibble(target_depth=c(30, 60), sd_log = c(.075, .075)))
  
  scenarios = scenarios0 %>%
    nest_by(scenario, target_depth) %>%
    mutate(bd_dist = list(discretize_bd_dist(slices, mean=target_depth+data$dd, data$sd_log))) %>% 
    unnest() %>%
    transmute(scenario, target_depth, sample_depth_max=depth, weight)
  
  # when compositing all, just put all weight on the mean change??
  scenarios_compAll = scenarios0 %>% 
    transmute(scenario, target_depth, sample_depth_max=target_depth+dd, weight=1) 
  
  bind_rows(
    scenarios %>% expand_grid(comp=c("none", "bottom")),
    scenarios_compAll %>% mutate(comp="all")
  )
}