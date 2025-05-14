# Load and install required packages, also need to install external rjags
library(MixSIAR)
library(here)
library(dplyr)
library(MixSIAR)
library(rjags)
library(MCMCpack)
library(GGally)
library(tidyverse)
library(ggplot2)
library(tidyr)


# Load mixture data (samples)
mix <- load_mix_data(
  filename = here("data","mixsiar_mixture_data.csv"),
  iso_names = c("d13C", "d15N"), #any isotopes run
  factors = "Site",       # add site as a fixed effect
  fac_random = FALSE,
  fac_nested = FALSE,
  cont_effects = NULL
)

# Load source/endmember data 
source <- load_source_data(
  filename = here("data","mixsiar_source_data.csv"),
  source_factors = NULL,
  conc_dep = FALSE,
  data_type = "means",
  mix = mix
)
# Load TEF (trophic enrichment factors) data (dummy vals for sediment) but need them otherwise program won't run
discr <- load_discr_data(
  filename = here("data","mixsiar_tef_data.csv"),
  mix = mix
)

# Plot tracer and sample data
plot_data(filename = "mixsiar_sediment_isotope_plot",
          mix = mix, source = source, discr = discr)

# Write the model
model_filename <- "mixsiar_model_sediment.txt"
write_JAGS_model(model_filename, resid_err = TRUE, process_err = TRUE, #resid is var in mixture, process is var in source
                 mix = mix, source = source)
#Run the model and wait a hella long time
jags.1 <- run_model(run = "normal",
                    mix = mix,
                    source = source,
                    discr = discr,
                    model_filename = model_filename)

#still waiting
##Observed stochastic nodes = the isotope data.
#Unobserved stochastic nodes = the unknown source proportions and other parameters.
#Total graph size = all parts the model is using to solve the mixing problem.

# Save the output results
output_options <- list(
  summary_save = TRUE,              
  summary_name = "mixsiar_summary", #posterior estimates
  sup_post = TRUE,                  
  plot_post_save_pdf = TRUE,        
  plot_post_name = "mixsiar_posterior_density", 
  plot_pairs = TRUE,                
  plot_pairs_save_pdf = TRUE,       
  plot_pairs_name = "mixsiar_pairs",
  diag_save = TRUE,                 
  diag_name = "mixsiar_diagnostics",
  indiv_plot = FALSE,               
  plot_save_png = FALSE,            
  summary_save_csv = TRUE,          
  sup_pairs = FALSE,                
  sup_post_indiv = FALSE,           
  sup_pairs_indiv = FALSE,        
  sup_xy = FALSE                    
)

#summarize model output
output_JAGS(jags.1, mix, source, output_options)
post_samples <- as.mcmc(jags.1$BUGSoutput$sims.matrix)

# Extract posterior samples from mcmc
post_samples <- as.mcmc(jags.1$BUGSoutput$sims.matrix)
post_df <- as.data.frame(post_samples)

# Get only proportions from model output and rename hella columns 
prop_df <- post_df %>% dplyr::select(starts_with("p"))

prop_long <- prop_df %>%
  pivot_longer(cols = everything(),
               names_to = "Source_Site",
               values_to = "Proportion") %>%
  mutate(
    Source = case_when(
      grepl("\\[.*,1\\]", Source_Site) ~ "Mangrove",
      grepl("\\[.*,2\\]", Source_Site) ~ "Seagrass",
      grepl("\\[.*,3\\]", Source_Site) ~ "Terrestrial",
      TRUE ~ "Unknown"
    ),
    Site = case_when(
      grepl("fac1\\[1,", Source_Site) ~ "Site_A",
      grepl("fac1\\[2,", Source_Site) ~ "Site_B",
      grepl("fac1\\[3,", Source_Site) ~ "Site_C",
      grepl("global\\[", Source_Site) ~ "Global",
      TRUE ~ "Unknown"
    )
  ) %>%
  dplyr::select(Site, Source, Proportion)


# Summarize mean proportions of samples by site and 95% CI
prop_summary <- prop_long %>%
  filter(Site %in% c("Site_A", "Site_B", "Site_C")) %>%
  group_by(Site, Source) %>%
  summarise(mean_prop = mean(Proportion),
            lower = quantile(Proportion, 0.025),
            upper = quantile(Proportion, 0.975),
            .groups = "drop")


# Plot stacked bar chart
stacked_source_cont<-ggplot(prop_summary, aes(x = Site, y = mean_prop, fill = Source)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Site", y = "Mean Source Contribution (%)") +
  theme_bw() +
  theme(legend.title = element_blank())
stacked_source_cont
#save
ggsave(here("plots", "source_cont.png"), plot=source_cont, width = 5, height=5, dpi=300)




##############################################################################################
#Part 2: run with SIMMR
install.packages("simmr")
library(simmr)
#load in mixture and source data
mixture_df <- read.csv(here("data","mixsiar_mixture_data_formatted.csv")) 
source_df <- read.csv(here("data", "mixsiar_source_data_corrected.csv"))

###plot c vs n like we did in mixsiar
ggplot() +
  geom_point(data = mixture_df, aes(x = d13C, y = d15N, color = Site), size = 3) +
  geom_point(data = source_df, aes(x = d13C, y = d15N), shape = 17, size = 5, color = "black") +
  geom_text(data = source_df, aes(x = d13C, y = d15N, label = Source), hjust = -0.2, vjust = -0.5) +
  theme_bw() +
  labs(title = "δ13C vs δ15N: Sediment Samples and Sources",
       x = expression(delta^13*C~("\u2030")),
       y = expression(delta^15*N~("\u2030")))


# Run simmr per site as for loop
for (site in unique(mixture_df$Site)) {
  site_data <- filter(mixture_df, Site == site)
  
  simmr_in <- simmr_load(
    mixtures = site_data[, c("d13C", "d15N")],
    source_names = source_df$Source,
    source_means = as.matrix(source_df[, c("Meand13C", "Meand15N")]),
    source_sds = as.matrix(source_df[, c("SDd13C", "SDd15N")]),
    correction_means = matrix(0, nrow = 3, ncol = 2),
    correction_sds = matrix(0, nrow = 3, ncol = 2)
  )
  
  simmr_out <- simmr_mcmc(simmr_in)
  results[[site]] <- simmr_out
}  
summary(results[["Site_B"]])  # or Site_B, Site_C
plot(results[["Site_B"]], type = "boxplot")


# Create a summary df of mean contributions for each source and site
summary_df <- lapply(names(results), function(site) {
  sim_summary <- summary(results[[site]], type = "statistics")$statistics
###plot as stacked bars
site_summary<-ggplot(summary_df, aes(x = Site, y = mean_prop, fill = Source)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Site", y = "Mean Source Contribution") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"))
site_summary

#load c:n data
c_n_ratio<-read.csv(here("data","c_n_ratio.csv"))
# Create the source_boxes data frame
source_boxes <- data.frame(
  Source = c("Terrestrial", "Seagrass", "Mangrove"),
  ymin = c(-27 - 1.0, -20 - 1.5, -24 - 1.0),
  ymax = c(-27 + 1.0, -20 + 1.5, -24 + 1.0),
  xmin = c(13.07 - 1.5, 10.94 - 1.5, 8.43 - 1.5),
  xmax = c(17.07 + 1.5, 10.94 + 1.5, 8.43 + 1.5)
)
source_boxes$label_x <- (source_boxes$xmin + source_boxes$xmax)/2
source_boxes$label_y <- source_boxes$ymax + 0.5
  

# Your plot with points and source boxes
c_n<-ggplot(c_n_ratio, aes(x = CN_ratio, y = d13C, color = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_rect(data = source_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", linetype = "dashed", inherit.aes = FALSE) +
  geom_text(data = source_boxes,
            aes(x = label_x, y = label_y, label = Source),
            inherit.aes = FALSE, size = 4)+
  labs(x = "C:N", 
       y = "d13C", 
       color = "Site") +
  theme_classic()
c_n


