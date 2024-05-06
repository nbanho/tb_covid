#### Libraries ####

library(tidyverse)


#### Raw data ####

raw_df <- read.csv("data-raw/TB_COVID19_survey_FUP_2023_cleaned.csv") %>%
  filter(country_ap != "Japan")

#### Site characteristics ####

site <- raw_df %>%
  mutate(
    site_region = factor(
      region,
      levels = c(
        "Asia-Pacific",
        "Central Africa",
        "East Africa",
        "Southern Africa",
        "West Africa"
      )
    )
  ) %>%
  mutate(
    site_country = trimws(paste0(
      country_ap,
      country_ca,
      country_ea,
      country_sa,
      country_wa
    )),
    site_country = factor(site_country, levels = sort(unique(site_country)))
  ) %>%
  mutate(
    site_setting = ifelse(grepl("Peri", clinic_setting),
      "peri-urban", clinic_setting
    ),
    site_setting = tolower(site_setting),
    site_setting = factor(
      site_setting,
      levels = c("urban", "peri-urban", "rural")
    )
  ) %>%
  mutate(
    site_owner = ifelse(grepl("Public", clinic_own), "public", "other"),
    site_owner = factor(site_owner, levels = c("public", "other"))
  ) %>%
  mutate(
    site_carelevel = stringi::stri_extract(level_of_care, regex = "\\w*"),
    site_carelevel = ifelse(site_carelevel == "Teritary",
      "Tertiary", site_carelevel
    ),
    site_carelevel = tolower(site_carelevel),
    site_carelevel = factor(
      site_carelevel,
      levels = c("primary", "secondary", "tertiary", "other")
    )
  ) %>%
  mutate(
    site_patients = ifelse(grepl("Both", patients), "in and out", "only out"),
    site_patients = factor(site_patients, levels = c("in and out", "only out"))
  ) %>%
  mutate(
    site_treated = ifelse(grepl("Both", treated),
      "adults and child", "only adults"
    ),
    site_treated = factor(
      site_treated,
      levels = c("adults and child", "only adults")
    )
  ) %>%
  mutate(
    site_integration =
      ifelse(grepl("Other", integration), "partial",
        ifelse(grepl("Full-integration", integration), "full",
          ifelse(grepl("Partial", integration), "partial", "not integrated")
        )
      ),
    site_integration = factor(
      site_integration,
      levels = c("full", "partial", "not integrated")
    )
  ) %>%
  dplyr::select(matches("site_"))

#### IPC measures ####

ipc <- raw_df %>%
  dplyr::select(
    matches("distance_"),
    matches("envmt_"),
    matches("msr_")
  ) %>%
  mutate_all(
    ~ ifelse(grepl("Already in place", .x), "before",
      ifelse(grepl("Only in place", .x), "during",
        ifelse(grepl("maintained", .x), "since",
          ifelse(grepl("Never", .x), "never",
            ifelse(grepl("Information", .x), "unknown", "not applicable")
          )
        )
      )
    )
  ) %>%
  mutate_all(
    ~ factor(.x, levels = c(
      "never", "before",
      "during", "since",
      "unknown", "not applicable"
    ))
  ) %>%
  dplyr::select(
    -distance_other,
    -distance_other_spec,
    -envmt_other,
    -envmt_other_spec,
    -msr_staff_other,
    -msr_staff_other_spec,
    -msr_patients_other,
    -msr_patients_other_spec
  )

#### Changes care provision #####

changes <- raw_df %>%
  dplyr::select(
    matches("changes_"),
    -matches("changes_dot"),
    -dot_changes_other_spec
  ) %>%
  mutate_all(
    ~ ifelse(grepl("Only during", .x), "during",
      ifelse(grepl("Started during", .x), "since",
        ifelse(grepl("Started after", .x), "after",
          ifelse(grepl("Not observed", .x), "never",
            ifelse(grepl("Information", .x), "unknown", "not applicable")
          )
        )
      )
    )
  ) %>%
  mutate_all(
    ~ factor(.x, levels = c(
      "never", "during", "since",
      "after", "unknown", "not applicable"
    ))
  ) %>%
  dplyr::select(
    -changes_other,
    -changes_other_spec,
    -changes_res_other,
    -matches("changes_res_spec")
  )

#### Changes DOT provision ####

changes_dot <- raw_df %>%
  dplyr::select(matches("changes_dot")) %>%
  mutate_all(
    ~ ifelse(grepl("Only in place", .x), "during",
      ifelse(grepl("Introduced during", .x), "since",
        ifelse(grepl("Never", .x), "never",
          ifelse(grepl("Information", .x), "unknown", "not applicable")
        )
      )
    )
  ) %>%
  mutate_all(
    ~ factor(.x, levels = c(
      "never", "during", "since",
      "unknown", "not applicable"
    ))
  ) %>%
  dplyr::select(-changes_dot_other)



#### Drug shortages ####

short <- raw_df %>%
  dplyr::select(
    matches("short_drugs_")
  ) %>%
  mutate_all(
    ~ ifelse(.x == "Yes", "yes", ifelse(.x == "No", "no", "not applicable"))
  ) %>%
  mutate_all(
    ~ factor(.x, levels = c("yes", "no", "not applicable"))
  )


#### Restrictions ####

restr <- raw_df %>%
  dplyr::select(
    matches("restr_"),
    -matches("tb_restr")
  ) %>%
  mutate_all(
    ~ ifelse(.x == "Yes", "yes", ifelse(.x == "No", "no", "not applicable"))
  ) %>%
  mutate_all(
    ~ factor(.x, levels = c("yes", "no", "not applicable"))
  ) %>%
  dplyr::select(
    -restr_tb_other
  )

#### Restrictions addressed ####

restr_addr <- raw_df %>%
  dplyr::select(matches("tb_restr_addr")) %>%
  mutate_all(
    ~ ifelse(grepl("Only in place", .x), "during",
      ifelse(grepl("Introduced during", .x), "since",
        ifelse(grepl("Never", .x), "never",
          ifelse(grepl("Information", .x), "unknown", "not applicable")
        )
      )
    )
  ) %>%
  mutate_all(
    ~ factor(.x, levels = c(
      "never", "during", "since",
      "unknown", "not applicable"
    ))
  ) %>%
  dplyr::select(-tb_restr_addr_other)


#### Combine characteristics ####

# combine
site_descr <- cbind(
  record_id = 1:nrow(raw_df),
  site, ipc, changes,
  changes_dot, short,
  restr, restr_addr
) %>%
  mutate(
    across(
      matches("restr_tb_"),
      ~ ifelse(site_integration == "not integrated",
        "not applicable", as.character(.x)
      )
    ),
    across(
      matches("restr_tb_"),
      ~ factor(.x, levels = c("yes", "no", "not applicable"))
    )
  )

# save
saveRDS(site_descr, "data-clean/site-characteristics.rds")

#### Case outcomes ####

# long format
df_long <- raw_df %>%
  dplyr::select(
    record_id,
    matches("_new_"),
    matches("_art_"),
    matches("_tpt_"),
    matches("_genex_")
  ) %>%
  mutate(record_id = 1:n()) %>%
  reshape2::melt("record_id") %>%
  mutate(
    quarter = as.numeric(
      gsub(
        "q", "",
        stringi::stri_extract(variable, regex = "q\\d")
      )
    ),
    year = as.numeric(
      paste0("20", gsub(
        "_", "",
        stringi::stri_extract(variable, regex = "_\\d{2}")
      ))
    ),
    outcome = gsub("_", "", stringi::stri_extract(variable, regex = "_\\w*_"))
  )

# quarterly data
df_quarterly <- df_long %>%
  dplyr::select(-variable) %>%
  mutate(year_quart = year + quarter / 4) %>%
  reshape2::dcast(record_id + year + quarter + year_quart ~ outcome) %>%
  mutate(across(c(art, genex, new, tpt), as.numeric))

# only annual data
df_annual <- df_quarterly %>%
  filter(is.na(quarter)) %>%
  dplyr::select(-quarter)

# filter years in quarterly data
df_quarterly <- df_quarterly %>%
  filter(!is.na(quarter))

# save data
saveRDS(df_quarterly, "data-clean/quarterly-data.rds")
saveRDS(df_annual, "data-clean/annual-data.rds")
