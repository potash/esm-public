# functions to help read in the finely sliced resample data

get_resample_labs = function() {

  lab_files = list.files("~/uiuc/esm/data/resample/", 
                         "TCTN",
                         full.names=TRUE)
  
  labs = do.call(bind_rows, lapply(lab_files, function(filename) { 
    read_csv(filename) %>% 
      mutate(filename=filename) %>%
      rename(any_of(c(`Sample ID`="...5"))) %>%
      rename(any_of(c(`Sample ID`="ID")))
  })) %>%
    mutate(`Sample ID` = str_replace(`Sample ID`, "OgleC", "Ogle5C")) %>%
    filter(!is.na(`Sample ID`)) %>%
    rename(`TOC%`=`C  [%]`) %>%
    filter(!(grepl("7.23.2025", filename) & `No.` >= 93)) %>%
    filter(!(grepl("6.18.2025", filename) & `No.` == 38))
  
  labs_costech = readxl::read_excel("~/uiuc/esm/data/resample/Eric_Lenarth_TC_TOC_TIC.xlsx", sheet=3) %>%
    select(-`Sample ID`) %>%
    rename(`Sample ID` = `...1`)
  
  r = bind_rows(labs, labs_costech) %>%
    extract(`Sample ID`, c("site_id", "core_id", "sample_depth_min", "sample_depth_max"),
            "LegP[_-](.*)[_-](.*)[_-](.*)-(.*)", remove=FALSE, convert=TRUE)
  
  # manually reran one sample to correct TIC%
  r %>%
    mutate(`TIC%` = ifelse(`Sample ID` == "LegP_Chri10_4_62.5-65", 4.2768,`TIC%`),
           `TOC%` = ifelse(`Sample ID` == "LegP_Chri10_4_62.5-65", `TC%` - `TIC%`, `TOC%`))
}

get_resample_soils = function() {
  soil_files = list.files("~/uiuc/esm/data/resample/", 
                        "Preliminary",
                        full.names=TRUE)

  soils = do.call(bind_rows, lapply(soil_files, function(filename) { 
    readxl::read_excel(filename) %>% 
      mutate(filename=filename) %>%
      mutate(`Total length (cm)` = as.numeric(`Total length (cm)`))
  })) %>%
    rename(site_id=`Field site:`,
           core_id = `Core:`) %>%
    extract(`Depth(cm):`, c("sample_depth_min", "sample_depth_max"),
            "(.*)-(.*)", remove=FALSE, convert=TRUE) %>%
    filter(!is.na(`Air dried soil (g)`))
  
  soils
}

get_resample_measurements = function(labs, soils) {
  soils_BD =  soils %>%
    mutate(thickness = sample_depth_max - sample_depth_min) %>%
    mutate(volume = thickness * pi * 2.2^2) %>%
    mutate(BD_air = `Air dried soil (g)`/volume) %>%
    mutate(Oven_to_air_ratio = `oven dried soil mass`/(`Air dried soil + tin (g)` -`Tin (g)`)) %>%
    mutate(Weight_oven = `Air dried soil (g)`*Oven_to_air_ratio) %>%
    mutate(BD = Weight_oven/volume)

  labs_soils = labs %>%
    select(-filename) %>%
    left_join(soils_BD %>% select(-`Sample ID`))

  labs_soils %>%
    mutate(location_id = paste0(site_id, "_", core_id), .after=site_id) %>%
    # this core is weird
    #filter((`site_id` == "Craw3" & core_id==11)) %>% View
    transmute(location_id, 
            sample_depth_min, sample_depth_max, 
            SOCc=`TOC%`,
            BD=BD) %>%
    mutate(location_id = str_replace(location_id, "Ogle5C", "Ogle5"))
}

