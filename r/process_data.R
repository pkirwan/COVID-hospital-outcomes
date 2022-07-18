# R file to process the extracted and linked data

here::i_am("r/data_processing.R")

library(tidyverse)
library(lubridate)
library(slider)
library(here)

date_data <- as_date("2021-11-22")
end_date <- as_date("2021-09-30")
level_order <- c(
    "Mar 2020", "Apr 2020", "May 2020", "Jun 2020", "Jul 2020", "Aug 2020", "Sep 2020",
    "Oct 2020", "Nov 2020", "Dec 2020", "Jan 2021", "Feb 2021", "Mar 2021", "Apr 2021",
    "May 2021", "Jun 2021", "Jul 2021", "Aug 2021", "Sep 2021"
)

sfp_dat <- readRDS(here("data/imported_data.rds"))

sfp_dat <- sfp_dat %>% filter(nosocomial == 0)

## bed occupancy data
bed_occupancy <- sfp_dat %>%
    filter(
        !is.na(hospital_in),
        hospital_in >= as_date("2020-03-01"),
        hospital_in <= end_date
    ) %>%
    mutate(
        date = case_when(
            specimen_date <= hospital_in ~ hospital_in,
            hospital_in < specimen_date ~ specimen_date
        )
    ) %>%
    arrange(date) %>%
    group_by(trust_name, date) %>%
    summarise(n = n()) %>%
    mutate(rollsum = slide_dbl(n, ~ sum(.x), .before = 3, .after = 3)) %>%
    mutate(pct = round(n / sum(n) * 100, 1))

bed_occupancy_rate <- bed_occupancy %>%
    group_by(trust_name) %>%
    slice_max(rollsum) %>%
    distinct(trust_name, .keep_all = TRUE) %>%
    select(trust_name, N = rollsum) %>%
    full_join(bed_occupancy, by = "trust_name") %>%
    mutate(bed_occupancy_rate = round(rollsum / N * 100, 1)) %>%
    select(trust_name, bed_occupancy_rate, date)

## Create the hosp datasets to use for estimates
## Assume those with unknown final outcome still alive, since data linked to death register

hosp <- sfp_dat %>%
    select(
        hospital_in, hospital_out, linkset, age, agegrp, sex, date_of_death, phec_name, utla_name,
        ethnicity_final, sgtf_under30CT, specimen_date, First, Second, Third, imd_rank,
        imd_decile, charlson_index, trust_name, provider_code, specimen_date, pillar
    ) %>%
    filter(
        !is.na(hospital_in),
        hospital_in >= as_date("2020-03-01"),
        hospital_in <= end_date
    ) %>%
    filter(
        hospital_out <= date_data | is.na(hospital_out),
        !is.na(phec_name)
    ) %>%
    mutate(
        date_of_death = dmy(date_of_death),
        hosp_interval = as.numeric(hospital_out - hospital_in),
        ttdeath = as.numeric(date_of_death - hospital_in),
        censor_date = as.numeric(date_data - hospital_in),
        censor_date = if_else(censor_date > 90, 90, censor_date),

        # definition for in-hospital death
        eventm = case_when(
            !is.na(date_of_death) & as.numeric(date_of_death - specimen_date) <= 90 & linkset == "CoV:ECDS" ~ "Death",
            !is.na(date_of_death) & as.numeric(date_of_death - specimen_date) <= 90 & linkset != "CoV:ECDS" & as.numeric(date_of_death - hospital_out) <= 14 ~ "Death",
            as.numeric(hospital_out - specimen_date) <= 90 & linkset != "CoV:ECDS" ~ "Discharge",
            TRUE ~ "Unknown final outcome"
        ),

        ## Event or right censoring time.
        timem = case_when(
            eventm == "Death" ~ ttdeath,
            eventm == "Discharge" ~ hosp_interval,
            !is.na(hosp_interval) & (hosp_interval < ttdeath | is.na(ttdeath)) & linkset != "CoV:ECDS" ~ hosp_interval,
            TRUE ~ 0.5
        ),

        ## Censoring status. Either exact+event known, interval censored + event known, or right censored  + event unknown
        ## Unknown outcomes considered right censored (used for times to death or ICU)
        statusm = if_else(eventm %in% c("Death", "Discharge"), 1, 3),

        ## Event time for time to death or time to ICU components.
        ## Exact time (time1,time1), right cens (time1,Inf), or int cens (0.5, time2)
        time1m = timem,
        time2m = if_else(statusm == 3, Inf, timem),
        time2m = if_else(time2m == Inf & !is.na(date_of_death), ttdeath, time2m),
        time1m = if_else(time1m == 0, 0.5, time1m),
        time1m = if_else(time1m > 90 & eventm == "Unknown final outcome", 90, time1m),
        time2m = if_else(time2m == 0, 0.5, time2m),
        eventm = factor(eventm, levels = c("Death", "Discharge"), ordered = FALSE),
        sex = factor(sex, levels = c("Male", "Female")),

        ## factor groups
        ageGrp7 = cut(age, c(0, 15, 25, 45, 65, 75, 85, Inf), right = FALSE),
        ageGrp7 = factor(ageGrp7, labels = c("0-14", "15-24", "25-44", "45-64", "65-74", "75-84", "85+")),
        pheRegion = factor(phec_name, levels = c("London", "East Midlands", "East of England", "North East", "North West", "South East", "South West", "West Midlands", "Yorkshire and Humber"), ordered = FALSE),
        ethnicity_final = factor(ethnicity_final, ordered = FALSE),
        ethGrp4 = fct_collapse(ethnicity_final,
            Asian = c("Any other Asian background", "Bangladeshi (Asian or Asian British)", "Indian (Asian or Asian British)", "Indian (Asian or Asian British)", "Pakistani (Asian or Asian British)"),
            Black = c("African (Black or Black British)", "Any other Black background", "Caribbean (Black or Black British)"),
            `Mixed/Other/Unknown` = c("Any other ethnic group", "Any other Mixed background", "Chinese (other ethnic group)", "White and Asian (Mixed)", "White and Black African (Mixed)", "White and Black Caribbean (Mixed)", "Unknown"),
            White = c("Any other White background", "British (White)", "Irish (White)")
        ),
        ethGrp4 = fct_relevel(ethGrp4, "White", "Asian", "Black", "Mixed/Other/Unknown"),
        monthyear = paste0(lubridate::month(hospital_in, label = TRUE), " ", year(hospital_in)),
        monthyear = factor(monthyear, levels = level_order),
        vaccine = case_when(
            specimen_date - 14 >= Second ~ "\u226514 days post second dose",
            specimen_date - 21 >= First ~ "\u226521 days post first dose",
            specimen_date >= First ~ "<21 days post First dose"
        ),
        vaccine = if_else(is.na(vaccine), " No vaccine", vaccine),
        imd_quintile = ceiling(imd_decile / 2) %>% as.character() %>% factor(),
        utla_name = factor(utla_name),
        charlson_index = factor(charlson_index, levels = c("0", "1-2", "3-4", ">=5"), labels = c("0", "1-2", "3-4", "5+")),
        trust_name = factor(trust_name, ordered = FALSE),
        date = case_when(
            specimen_date <= hospital_in ~ hospital_in,
            hospital_in < specimen_date ~ specimen_date
        ),
        week = paste0(year(specimen_date), "-", if_else(lubridate::week(specimen_date) < 10, "0", ""), lubridate::week(specimen_date))
    ) %>%
    filter(!is.na(time1m), time1m >= 0, time1m <= 90) %>%
    # join the bed_occupancy rate data
    left_join(bed_occupancy_rate, by = c("trust_name", "date")) %>%
    mutate(
        bed_occupancy_cat = cut(bed_occupancy_rate, breaks = c(0, 20, 40, 60, 80, 90, Inf), right = FALSE),
        bed_occupancy_cat = factor(bed_occupancy_cat,
            levels = c("[0,20)", "[20,40)", "[40,60)", "[60,80)", "[80,90)", "[90,Inf)"),
            labels = c("0-20", "20-40", "40-60", "60-80", "80-90", "90-100")
        ),
    )

save(hosp, level_order, date_data, end_date, file = here("data/processed_data.RData"))
