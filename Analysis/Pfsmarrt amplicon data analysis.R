
######################################################################
######################################################################
######################################################################
# Pfsmarrt data amplicon data drug resistance analysis
# By Abebe A. Fola
# Date 04042025 
######################################################################
######################################################################
######################################################################


#Start by telling your computer where the files/set the working directory
getwd()
setwd("C:/Users/afola/Desktop/EPHI_projects/outbreaksamples/newnextseq")
#Install Library

# The duration of microscopy positive over time 
# [time-to-becoming negative] was analyzed using the Cox-PHZ;
# with time of follow-up treated as the main covariate for evaluating the proportion of positive samples (by microscopy) as the outcome, with patient ID incorporated as a random effect [frailty function].

library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(tidyverse)
library(dplyr)
library(tidyr)
library(viridis)
library(leaflet)
library(rhandsontable)
library(sp)
library(rgeos)
library(ggpubr)
library(scales)
library(gridExtra)
library(scatterplot3d)
library(knitr)
library(ggExtra)
library(GGally)
library(viridis)
library("ggsci")

# Load necessary libraries
library(ggplot2)
library(RColorBrewer)


#Load AA fract data from seekdeep output 

AA_tabale <- read.cv ("amino_acid_fracs.csv") # select only drug resistance loci and add metadata as needed 
# Save as or write  DRgenotypes.csv and use the data for down stream analysis


#merge data drug, and metadata


data1<- read.csv("Outbreak_Metadata_Sequencedidcorrect.csv")

data2 <- read.csv("DRgenotypes.csv")

metadrmergedata<- left_join(data1, data2, by="Sample_ID")

write.csv(metadrmergedata, "metadrmergedata.csv")

metadrmergedata <- ("metadrmergedata.csv")

# Read in the data
data <- read.csv("metadrmergedata.csv")

# Check the column names in your data

# Read the CSV data (assuming it's saved as 'DRgenotypes.csv')
data <- read.csv("DRgenotypes.csv")

# Check the column names in your data
print(colnames(data))  # Check the exact column names in your dataset

# Define the markers vector (make sure it matches the exact column names in your data)
markers <- c("crt_76T", "dhfr_51I", "dhfr_59R", "dhfr_108N", "dhps_437G", "dhps_540E", "dhps_581G", 
             "k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V", 
             "mdr1_N86", "mdr1_184F", "mdr1_D1246")

# Define the function to calculate prevalence and confidence intervals
calculate_prevalence_with_ci <- function(marker_column) {
  # Remove NA values from the marker column
  cleaned_column <- marker_column[!is.na(marker_column)]
  
  # Classify as mutant (>0) or wildtype (0)
  marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
  
  # Calculate prevalence of 'mutant' status
  prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
  
  # Calculate the number of mutant and wildtype samples
  num_mutant <- sum(marker_classification == "mutant")
  num_wildtype <- sum(marker_classification == "wildtype")
  
  # Calculate 95% confidence intervals for prevalence (binomial CI)
  ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
  
  return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
           LowerCI = ci[1] * 100,  # multiply CI by 100
           UpperCI = ci[2] * 100,  # multiply CI by 100
           NumMutant = num_mutant,
           NumWildtype = num_wildtype))
}

# Initialize an empty list to store results
results_list <- list()

# Loop through each marker and calculate prevalence
for (marker in markers) {
  # Calculate prevalence and CI for the marker across all data
  marker_column <- data[[marker]]
  marker_results <- calculate_prevalence_with_ci(marker_column)
  
  # Store the results in a data frame
  results_df <- data.frame(Marker = marker,
                           Prevalence = marker_results["Prevalence"],
                           LowerCI = marker_results["LowerCI"],
                           UpperCI = marker_results["UpperCI"],
                           NumMutant = marker_results["NumMutant"],
                           NumWildtype = marker_results["NumWildtype"])
  
  # Add the results for this marker to the list
  results_list[[marker]] <- results_df
}

# Combine all marker results into one final data frame
results_df_final <- do.call(rbind, results_list)

# Write the results to a CSV file
write.csv(results_df_final, "prevalence_and_ci_across_all_sites.csv", row.names = FALSE)

# Create the bar plot with different colors for each marker (Color Blind-Friendly Palette)
ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Marker)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
  scale_fill_brewer(palette = "Set1") +  # Use colorblind-friendly colors
  theme_minimal() +
  labs(title = "Prevalence of Mutants for Each Marker with 95% CI",
       x = "Markers",
       y = "Prevalence of Mutants (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability


###
# gene level color
###

# Load necessary libraries
library(ggplot2)
library(RColorBrewer)

# Read the CSV data (assuming it's saved as 'DRgenotypes.csv')
data1 <- read.csv("DRgenotypes.csv")

# Check the column names in your data
print(colnames(data1))  # Check the exact column names in your dataset

# Define the markers vector (make sure it matches the exact column names in your data)
markers <- c("crt_76T", "dhfr_51I", "dhfr_59R", "dhfr_108N", "dhps_437G", "dhps_540E", "dhps_581G", 
             "k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V", 
             "mdr1_N86", "mdr1_184F", "mdr1_D1246")

# Create a custom color mapping based on marker families (initials)
marker_colors <- c(
  "crt" = "#e6194b",  # red for crt family
  "dhfr_ts" = "#f58231",  # orange for dhfr family
  "dhps" =  "#56B4E9" , # purple for dhps family
  "mdr1" = "#009E73",  # purple for mdr1 family
  "k13" =  "#CC79A7" # cyan for k13 family
)

# Map each marker to its family group (e.g., crt, dhfr_ts, dhps, mdr1, k13)
marker_family_map <- c(
  "crt_76T" = "crt",
  "dhfr_51I" = "dhfr_ts", "dhfr_59R" = "dhfr_ts", "dhfr_108N" = "dhfr_ts",
  "dhps_437G" = "dhps", "dhps_540E" = "dhps", "dhps_581G" = "dhps",
  "k13_433D" = "k13", "k13_441L" = "k13", "k13_574L" = "k13", "k13_622I" = "k13", 
  "k13_658E" = "k13", "k13_675V" = "k13",
  "mdr1_N86" = "mdr1", "mdr1_184F" = "mdr1", "mdr1_D1246" = "mdr1"
)

# Define the function to calculate prevalence and confidence intervals
calculate_prevalence_with_ci <- function(marker_column) {
  # Remove NA values from the marker column
  cleaned_column <- marker_column[!is.na(marker_column)]
  
  # Classify as mutant (>0) or wildtype (0)
  marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
  
  # Calculate prevalence of 'mutant' status
  prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
  
  # Calculate the number of mutant and wildtype samples
  num_mutant <- sum(marker_classification == "mutant")
  num_wildtype <- sum(marker_classification == "wildtype")
  
  # Calculate 95% confidence intervals for prevalence (binomial CI)
  ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
  
  return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
           LowerCI = ci[1] * 100,  # multiply CI by 100
           UpperCI = ci[2] * 100,  # multiply CI by 100
           NumMutant = num_mutant,
           NumWildtype = num_wildtype))
}

# Initialize an empty list to store results
results_list <- list()

# Loop through each marker and calculate prevalence
for (marker in markers) {
  # Calculate prevalence and CI for the marker across all data
  marker_column <- data1[[marker]]
  marker_results <- calculate_prevalence_with_ci(marker_column)
  
  # Store the results in a data frame
  results_df <- data.frame(Marker = marker,
                           Family = marker_family_map[marker],  # Assign family based on the map
                           Prevalence = marker_results["Prevalence"],
                           LowerCI = marker_results["LowerCI"],
                           UpperCI = marker_results["UpperCI"],
                           NumMutant = marker_results["NumMutant"],
                           NumWildtype = marker_results["NumWildtype"])
  
  # Add the results for this marker to the list
  results_list[[marker]] <- results_df
}

# Combine all marker results into one final data frame
results_df_final <- do.call(rbind, results_list)

# Write the results to a CSV file
write.csv(results_df_final, "prevalence_and_ci_by_family.csv", row.names = FALSE)

# Updated ggplot code to use "Genes" as the legend title
ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Family)) +
  geom_bar(stat = "identity") +  # Plot the bars
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
  scale_fill_manual(values = marker_colors) +  # Apply custom colors
  theme_minimal() +
  labs(title = "Overall Prevalence of Drug Resistance Markers",
       x = "Markers",
       y = "Prevalence of Mutants (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(legend.title = element_text(size = 12),          # Customize legend title font size
        legend.text = element_text(size = 10)) +          # Customize legend text font size
  guides(fill = guide_legend(title = "Genes"))  # Change the legend title to "Genes"



# Read the CSV data (assuming it's saved as 'DRgenotypes.csv') subset samples per area or study type
phemsamplesonly <- read.csv("phemsamplesonly.csv")

# Check the column names in your data
print(colnames(phemsamplesonly))  # Check the exact column names in your dataset

# Define the markers vector (make sure it matches the exact column names in your data)
markers <- c("crt_76T", "dhfr_51I", "dhfr_59R", "dhfr_108N", "dhps_437G", "dhps_540E", "dhps_581G", 
             "k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V", 
             "mdr1_N86", "mdr1_184F", "mdr1_D1246")

# Create a custom color mapping based on marker families (initials)
marker_colors <- c(
  "crt" = "#e6194b",  # red for crt family
  "dhfr_ts" = "#f58231",  # orange for dhfr family
  "dhps" =  "#56B4E9" , # purple for dhps family
  "mdr1" = "#009E73",  # purple for mdr1 family
  "k13" =  "#CC79A7" # cyan for k13 family
)

# Map each marker to its family group (e.g., crt, dhfr_ts, dhps, mdr1, k13)
marker_family_map <- c(
  "crt_76T" = "crt",
  "dhfr_51I" = "dhfr_ts", "dhfr_59R" = "dhfr_ts", "dhfr_108N" = "dhfr_ts",
  "dhps_437G" = "dhps", "dhps_540E" = "dhps", "dhps_581G" = "dhps",
  "k13_433D" = "k13", "k13_441L" = "k13", "k13_574L" = "k13", "k13_622I" = "k13", 
  "k13_658E" = "k13", "k13_675V" = "k13",
  "mdr1_N86" = "mdr1", "mdr1_184F" = "mdr1", "mdr1_D1246" = "mdr1"
)

# Define the function to calculate prevalence and confidence intervals
calculate_prevalence_with_ci <- function(marker_column) {
  # Remove NA values from the marker column
  cleaned_column <- marker_column[!is.na(marker_column)]
  
  # Classify as mutant (>0) or wildtype (0)
  marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
  
  # Calculate prevalence of 'mutant' status
  prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
  
  # Calculate the number of mutant and wildtype samples
  num_mutant <- sum(marker_classification == "mutant")
  num_wildtype <- sum(marker_classification == "wildtype")
  
  # Calculate 95% confidence intervals for prevalence (binomial CI)
  ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
  
  return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
           LowerCI = ci[1] * 100,  # multiply CI by 100
           UpperCI = ci[2] * 100,  # multiply CI by 100
           NumMutant = num_mutant,
           NumWildtype = num_wildtype))
}

# Initialize an empty list to store results
results_list <- list()

# Loop through each marker and calculate prevalence
for (marker in markers) {
  # Calculate prevalence and CI for the marker across all data
  marker_column <- phemsamplesonly[[marker]]
  marker_results <- calculate_prevalence_with_ci(marker_column)
  
  # Store the results in a data frame
  results_df <- data.frame(Marker = marker,
                           Family = marker_family_map[marker],  # Assign family based on the map
                           Prevalence = marker_results["Prevalence"],
                           LowerCI = marker_results["LowerCI"],
                           UpperCI = marker_results["UpperCI"],
                           NumMutant = marker_results["NumMutant"],
                           NumWildtype = marker_results["NumWildtype"])
  
  # Add the results for this marker to the list
  results_list[[marker]] <- results_df
}

# Combine all marker results into one final data frame
results_df_final <- do.call(rbind, results_list)

# Write the results to a CSV file
write.csv(results_df_final, "prevalence_and_ci_phemsamplesonly.csv", row.names = FALSE)

# Create the bar plot with custom color mapping based on marker families

# Updated ggplot code to use "Genes" as the legend title
ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Family)) +
  geom_bar(stat = "identity") +  # Plot the bars
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
  scale_fill_manual(values = marker_colors) +  # Apply custom colors
  theme_minimal() +
  labs(title = "Phemsamples-Prevalence of Drug Resistance Markers",
       x = "Markers",
       y = "Prevalence of Mutants (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(legend.title = element_text(size = 12),          # Customize legend title font size
        legend.text = element_text(size = 10)) +          # Customize legend text font size
  guides(fill = guide_legend(title = "Genes"))  # Change the legend title to "Genes"

#
# Work on site level
#

# Load necessary libraries
library(ggplot2)
getwd()
setwd("C:/Users/afola/Desktop/EPHI_projects/outbreaksamples/newnextseq")
# Read the CSV data (assuming it's saved as 'DRgenotypes.csv')
data1 <- read.csv("metadrmergedatacoimod.csv")
Bokre.Site = c("#e6194b",  "#f58231", "#ffe119",   "#911eb4", 
                        "#f032e6", "#4363d8","#46f0f0","#008080", "#3cb44b") 
                        barplotSite <- ggplot(data1, aes(site))
                        
                        barplotSite  + 
                          geom_bar(aes(fill=Site_Name)) +
                          #scale_fill_brewer(palette="set2")+
                          labs(x="Site", y="Number of samples") +
                          guides(fill=guide_legend(title="Sites"))+
                          #scale_fill_manual(values=Bokre.Site)+
                          #scale_fill_viridis_d()
                          theme(axis.text.x = element_text(angle = 90))+
                          theme() +
                          ggtitle(" Fig 2.Samples Distribution per Region per site") 
                        
                        
                        
                        # Load necessary library
                        library(dplyr)
                        
                        # Read the CSV file
                        df <- read.csv("Outbreak_Metadata_Sequencedidcorrect.csv")  # Replace with your actual file name
                        
                        # Create the summary table and order it
                        summary_table <- df %>%
                          group_by(Region, Woreda, Site_Name) %>%
                          summarise(Sample_Count = n(), .groups = "drop") %>%
                          arrange(Region, Woreda, Site_Name)
                        
                        # Save the summary to a new CSV
                        write.csv(summary_table, "summary_table.csv", row.names = FALSE)
                        
                        cat("Ordered summary table saved as 'summary_table.csv'\n")
                        
                        # Check the column names in your data
                        print(colnames(data1))  # Check the exact column names in your dataset
                        
                        # Define the markers vector (make sure it matches the exact column names in your data)
                        markers <- c("crt_76T", "dhfr_51I", "dhfr_59R", "dhfr_108N", "dhps_437G", "dhps_540E", "dhps_581G", 
                                     "k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V", 
                                     "mdr1_N86", "mdr1_184F", "mdr1_D1246")
                        
                        # Create a custom color mapping based on marker families (initials)
                        marker_colors <- c(
                          "crt" = "#e6194b",  # red for crt family
                          "dhfr_ts" = "#f58231",  # orange for dhfr family
                          "dhps" =  "#56B4E9", # purple for dhps family
                          "mdr1" = "#009E73",  # green for mdr1 family
                          "k13" =  "#CC79A7"  # cyan for k13 family
                        )
                        
                        # Map each marker to its family group (e.g., crt, dhfr_ts, dhps, mdr1, k13)
                        marker_family_map <- c(
                          "crt_76T" = "crt",
                          "dhfr_51I" = "dhfr_ts", "dhfr_59R" = "dhfr_ts", "dhfr_108N" = "dhfr_ts",
                          "dhps_437G" = "dhps", "dhps_540E" = "dhps", "dhps_581G" = "dhps",
                          "k13_433D" = "k13", "k13_441L" = "k13", "k13_574L" = "k13", "k13_622I" = "k13", 
                          "k13_658E" = "k13", "k13_675V" = "k13",
                          "mdr1_N86" = "mdr1", "mdr1_184F" = "mdr1", "mdr1_D1246" = "mdr1"
                        )
                        
                        # Define the function to calculate prevalence and confidence intervals
                        calculate_prevalence_with_ci <- function(marker_column) {
                          # Remove NA values from the marker column
                          cleaned_column <- marker_column[!is.na(marker_column)]
                          
                          # Classify as mutant (>0) or wildtype (0)
                          marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
                          
                          # Calculate prevalence of 'mutant' status
                          prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
                          
                          # Calculate the number of mutant and wildtype samples
                          num_mutant <- sum(marker_classification == "mutant")
                          num_wildtype <- sum(marker_classification == "wildtype")
                          
                          # Check if there are mutants to avoid error in binom.test
                          if (num_mutant > 0 && length(cleaned_column) > num_mutant) {
                            # If there are mutants and sufficient samples, run binom.test
                            ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
                            return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
                                     LowerCI = ci[1] * 100,  # multiply CI by 100
                                     UpperCI = ci[2] * 100,  # multiply CI by 100
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          } else {
                            # If no mutants or insufficient data, return NA for confidence intervals
                            return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
                                     LowerCI = NA,
                                     UpperCI = NA,
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          }
                        }
                        
                        # Initialize an empty list to store results
                        results_list <- list()
                        
                        # Loop through each site in your data
                        for (site in unique(data1$Site)) {
                          # Filter data for the current site
                          site_data <- data1[data1$Site == site, ]
                          
                          # Loop through each marker and calculate prevalence for each site
                          for (marker in markers) {
                            # Calculate prevalence and CI for the marker at the current site
                            marker_column <- site_data[[marker]]
                            marker_results <- calculate_prevalence_with_ci(marker_column)
                            
                            # Store the results in a data frame
                            results_df <- data.frame(Site = site,
                                                     Marker = marker,
                                                     Family = marker_family_map[marker],  # Assign family based on the map
                                                     Prevalence = marker_results["Prevalence"],
                                                     LowerCI = marker_results["LowerCI"],
                                                     UpperCI = marker_results["UpperCI"],
                                                     NumMutant = marker_results["NumMutant"],
                                                     NumWildtype = marker_results["NumWildtype"])
                            
                            # Add the results for this marker and site to the list
                            results_list[[paste(site, marker, sep = "_")]] <- results_df
                          }
                        }
                        
                        # Combine all marker results into one final data frame
                        results_df_final <- do.call(rbind, results_list)
                        
                        # Write the results to a CSV file
                        write.csv(results_df_final, "prevalence_and_ci_by_site_and_marker.csv", row.names = FALSE)
                        
                        # Create the bar plot with custom color mapping based on marker families
                        ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Family)) +
                          geom_bar(stat = "identity", show.legend = FALSE) +
                          geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
                          facet_wrap(~Site) +  # Add facets by site
                          scale_fill_manual(values = marker_colors) +  # Apply custom color mapping
                          theme_minimal() +
                          labs(title = "Prevalence of Mutants for Each Marker with 95% CI by Site",
                               x = "Markers",
                               y = "Prevalence of Mutants (%)") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
                        
                        # Updated ggplot code to use "Genes" as the legend title
                        ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Family)) +
                          geom_bar(stat = "identity") +  # Plot the bars
                          facet_wrap(~Site) +  # Add facets by site
                          geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
                          scale_fill_manual(values = marker_colors) +  # Apply custom colors
                          theme_minimal() +
                          labs(title = "Overall Prevalence of Drug Resistance Markers",
                               x = "Markers",
                               y = "Prevalence of Mutants (%)") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
                          theme(legend.title = element_text(size = 12),          # Customize legend title font size
                                legend.text = element_text(size = 10)) +          # Customize legend text font size
                          guides(fill = guide_legend(title = "Genes"))  # Change the legend title to "Genes"
                        
                        
                        # Regional level
                        
                        # Check the column names in your data
                        print(colnames(data1))  # Check the exact column names in your dataset
                        
                        # Define the markers vector (make sure it matches the exact column names in your data)
                        markers <- c("crt_76T", "dhfr_51I", "dhfr_59R", "dhfr_108N", "dhps_437G", "dhps_540E", "dhps_581G", 
                                     "k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V", 
                                     "mdr1_N86", "mdr1_184F", "mdr1_D1246")
                        
                        # Create a custom color mapping based on marker families (initials)
                        marker_colors <- c(
                          "crt" = "#e6194b",  # red for crt family
                          "dhfr_ts" = "#f58231",  # orange for dhfr family
                          "dhps" =  "#56B4E9", # purple for dhps family
                          "mdr1" = "#009E73",  # green for mdr1 family
                          "k13" =  "#CC79A7"  # cyan for k13 family
                        )
                        
                        # Map each marker to its family group (e.g., crt, dhfr_ts, dhps, mdr1, k13)
                        marker_family_map <- c(
                          "crt_76T" = "crt",
                          "dhfr_51I" = "dhfr_ts", "dhfr_59R" = "dhfr_ts", "dhfr_108N" = "dhfr_ts",
                          "dhps_437G" = "dhps", "dhps_540E" = "dhps", "dhps_581G" = "dhps",
                          "k13_433D" = "k13", "k13_441L" = "k13", "k13_574L" = "k13", "k13_622I" = "k13", 
                          "k13_658E" = "k13", "k13_675V" = "k13",
                          "mdr1_N86" = "mdr1", "mdr1_184F" = "mdr1", "mdr1_D1246" = "mdr1"
                        )
                        
                        # Define the function to calculate prevalence and confidence intervals
                        calculate_prevalence_with_ci <- function(marker_column) {
                          # Remove NA values from the marker column
                          cleaned_column <- marker_column[!is.na(marker_column)]
                          
                          # Classify as mutant (>0) or wildtype (0)
                          marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
                          
                          # Calculate prevalence of 'mutant' status
                          prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
                          
                          # Calculate the number of mutant and wildtype samples
                          num_mutant <- sum(marker_classification == "mutant")
                          num_wildtype <- sum(marker_classification == "wildtype")
                          
                          # Check if there are mutants to avoid error in binom.test
                          if (num_mutant > 0 && length(cleaned_column) > num_mutant) {
                            # If there are mutants and sufficient samples, run binom.test
                            ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
                            return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
                                     LowerCI = ci[1] * 100,  # multiply CI by 100
                                     UpperCI = ci[2] * 100,  # multiply CI by 100
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          } else {
                            # If no mutants or insufficient data, return NA for confidence intervals
                            return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
                                     LowerCI = NA,
                                     UpperCI = NA,
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          }
                        }
                        
                        # Initialize an empty list to store results
                        results_list <- list()
                        
                        # Loop through each site in your data
                        for (site in unique(data1$site)) {
                          # Filter data for the current site
                          site_data <- data1[data1$site == site, ]
                          
                          # Loop through each marker and calculate prevalence for each site
                          for (marker in markers) {
                            # Calculate prevalence and CI for the marker at the current site
                            marker_column <- site_data[[marker]]
                            marker_results <- calculate_prevalence_with_ci(marker_column)
                            
                            # Store the results in a data frame
                            results_df <- data.frame(site = site,
                                                     Marker = marker,
                                                     Family = marker_family_map[marker],  # Assign family based on the map
                                                     Prevalence = marker_results["Prevalence"],
                                                     LowerCI = marker_results["LowerCI"],
                                                     UpperCI = marker_results["UpperCI"],
                                                     NumMutant = marker_results["NumMutant"],
                                                     NumWildtype = marker_results["NumWildtype"])
                            
                            # Add the results for this marker and site to the list
                            results_list[[paste(site, marker, sep = "_")]] <- results_df
                          }
                        }
                        
                        # Combine all marker results into one final data frame
                        results_df_final <- do.call(rbind, results_list)
                        
                        # Write the results to a CSV file
                        write.csv(results_df_final, "prevalence_and_ci_by_site_and_marker.csv", row.names = FALSE)
                        
                        # Create the bar plot with custom color mapping based on marker families
                        # Updated ggplot code to use "Genes" as the legend title
                        ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Family)) +
                          geom_bar(stat = "identity") +  # Plot the bars
                          facet_wrap(~site) +  # Add facets by site
                          geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
                          scale_fill_manual(values = marker_colors) +  # Apply custom colors
                          theme_minimal() +
                          labs(title = "Overall Prevalence of Drug Resistance Markers",
                               x = "Markers",
                               y = "Prevalence of Mutants (%)") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
                          theme(legend.title = element_text(size = 12),          # Customize legend title font size
                                legend.text = element_text(size = 10)) +          # Customize legend text font size
                          guides(fill = guide_legend(title = "Genes"))  # Change the legend title to "Genes"
                        
                        
                        
                        #Only K13 mutations 
                        
                        # Load RColorBrewer for the Set1 palette
                        setwd("C:/Users/afola/Desktop/EPHI_projects/outbreaksamples/newnextseq")
                        
                        
                        library(RColorBrewer)
                        
                        data1 <- read.csv("metadrmergedata.csv")
                        
                        
                        # Define the k13 family markers
                        markers <- c("k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V")
                        
                        # Generate distinct colors using the Set1 palette
                        marker_colors <- brewer.pal(length(markers), "Set1")
                        
                        # Map each marker to its assigned color
                        marker_color_map <- setNames(marker_colors, markers)
                        
                        # Define the function to calculate prevalence and confidence intervals
                        calculate_prevalence_with_ci <- function(marker_column) {
                          # Remove NA values from the marker column
                          cleaned_column <- marker_column[!is.na(marker_column)]
                          
                          # Classify as mutant (>0) or wildtype (0)
                          marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
                          
                          # Calculate prevalence of 'mutant' status
                          prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
                          
                          # Calculate the number of mutant and wildtype samples
                          num_mutant <- sum(marker_classification == "mutant")
                          num_wildtype <- sum(marker_classification == "wildtype")
                          
                          # Check if there are mutants to avoid error in binom.test
                          if (num_mutant > 0 && length(cleaned_column) > num_mutant) {
                            # If there are mutants and sufficient samples, run binom.test
                            ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
                            return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
                                     LowerCI = ci[1] * 100,  # multiply CI by 100
                                     UpperCI = ci[2] * 100,  # multiply CI by 100
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          } else {
                            # If no mutants or insufficient data, return NA for confidence intervals
                            return(c(Prevalence = prevalence * 100,  # multiply by 100 to convert to percentage
                                     LowerCI = NA,
                                     UpperCI = NA,
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          }
                        }
                        
                        # Initialize an empty list to store results
                        results_list <- list()
                        
                        # Loop through each site in your data (only for the k13 markers)
                        for (site in unique(data1$site)) {
                          # Filter data for the current site
                          site_data <- data1[data1$site == site, ]
                          
                          # Loop through each k13 marker and calculate prevalence for each site
                          for (marker in markers) {
                            # Calculate prevalence and CI for the marker at the current site
                            marker_column <- site_data[[marker]]
                            marker_results <- calculate_prevalence_with_ci(marker_column)
                            
                            # Store the results in a data frame
                            results_df <- data.frame(site = site,
                                                     Marker = marker,
                                                     Family = "k13",  # All are k13 family markers
                                                     Prevalence = marker_results["Prevalence"],
                                                     LowerCI = marker_results["LowerCI"],
                                                     UpperCI = marker_results["UpperCI"],
                                                     NumMutant = marker_results["NumMutant"],
                                                     NumWildtype = marker_results["NumWildtype"])
                            
                            # Add the results for this marker and site to the list
                            results_list[[paste(site, marker, sep = "_")]] <- results_df
                          }
                        }
                        
                        # Combine all k13 marker results into one final data frame
                        results_df_final <- do.call(rbind, results_list)
                        
                        # Write the results to a CSV file
                        write.csv(results_df_final, "k13_prevalence_and_ci_by_site.csv", row.names = FALSE)
                        
                        # Create the bar plot with custom color mapping based on marker families
                        # Use distinct colors from Set1 for each marker
                        ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Marker)) +
                          geom_bar(stat = "identity") +  # Plot the bars
                          #facet_wrap(~site) +  # Add facets by site
                          geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
                          scale_fill_manual(values = marker_color_map) +  # Apply custom colors from Set1
                          theme_minimal() +
                          labs(title = "Outbreak Samples–Spatial Heterogeneity of K13 Mutations (ART-R Markers))",
                               x = "Markers",
                               y = "Prevalence of Mutants (%)") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
                          theme(legend.title = element_text(size = 12),          # Customize legend title font size
                                legend.text = element_text(size = 10)) +          # Customize legend text font size
                          guides(fill = guide_legend(title = "Markers"))  # Change the legend title to "Markers"
                        
                        
                        
                        
                        #####
                        # Over all k13 
                        ######
                        
                        # Load necessary libraries
                        library(dplyr)
                        library(ggplot2)
                        library(RColorBrewer)
                        
                        # Read your data
                        data1 <- read.csv("metadrmergedata.csv")
                        
                        
                        # Define the K13 markers
                        markers <- c("k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V")
                        
                        # Generate distinct colors using the Set1 palette
                        marker_colors <- brewer.pal(length(markers), "Set1")
                        marker_color_map <- setNames(marker_colors, markers)
                        
                        # Define function to calculate prevalence and confidence intervals
                        calculate_prevalence_with_ci <- function(marker_column) {
                          cleaned_column <- marker_column[!is.na(marker_column)]
                          marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
                          
                          prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
                          num_mutant <- sum(marker_classification == "mutant")
                          num_wildtype <- sum(marker_classification == "wildtype")
                          
                          if (num_mutant > 0 && length(cleaned_column) > num_mutant) {
                            ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
                            return(c(Prevalence = prevalence * 100,
                                     LowerCI = ci[1] * 100,
                                     UpperCI = ci[2] * 100,
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          } else {
                            return(c(Prevalence = prevalence * 100,
                                     LowerCI = NA,
                                     UpperCI = NA,
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          }
                        }
                        
                        # Initialize results list for overall prevalence
                        results_list <- list()
                        
                        # Loop through each marker (no site grouping)
                        for (marker in markers) {
                          marker_column <- data1[[marker]]
                          marker_results <- calculate_prevalence_with_ci(marker_column)
                          
                          results_df <- data.frame(Marker = marker,
                                                   Family = "k13",
                                                   Prevalence = marker_results["Prevalence"],
                                                   LowerCI = marker_results["LowerCI"],
                                                   UpperCI = marker_results["UpperCI"],
                                                   NumMutant = marker_results["NumMutant"],
                                                   NumWildtype = marker_results["NumWildtype"])
                          
                          results_list[[marker]] <- results_df
                        }
                        
                        # Combine results
                        results_df_final <- do.call(rbind, results_list)
                        
                        # Save to CSV
                        write.csv(results_df_final, "k13_prevalence_overall.csv", row.names = FALSE)
                        
                        # Plot bar 
                        ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Marker)) +
                          geom_bar(stat = "identity") +
                          geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
                          scale_fill_manual(values = marker_color_map) +
                          theme_minimal() +
                          labs(title = "K13 Mutation Prevalence (Overall)",
                               x = "Markers",
                               y = "Prevalence of Mutants (%)") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                legend.title = element_text(size = 12),
                                legend.text = element_text(size = 10)) +
                          guides(fill = guide_legend(title = "Markers"))
                        
                        
                        
                        
                        
                        ##########
                        # Site level k13 prevalence
                        #########
                        
                        library(ggplot2)
                        library(RColorBrewer)
                        library(dplyr)
                        
                        # Load the data
                        data1 <- read.csv("metadrmergedata.csv")
                        
                        # Define the K13 markers
                        markers <- c("k13_433D", "k13_441L", "k13_574L", "k13_622I", "k13_658E", "k13_675V")
                        
                        # Generate distinct colors using the Set1 palette
                        marker_colors <- brewer.pal(length(markers), "Set1")
                        marker_color_map <- setNames(marker_colors, markers)
                        
                        # Define function to calculate prevalence and confidence intervals
                        calculate_prevalence_with_ci <- function(marker_column) {
                          cleaned_column <- marker_column[!is.na(marker_column)]
                          marker_classification <- ifelse(cleaned_column > 0, "mutant", "wildtype")
                          
                          prevalence <- sum(marker_classification == "mutant") / length(cleaned_column)
                          num_mutant <- sum(marker_classification == "mutant")
                          num_wildtype <- sum(marker_classification == "wildtype")
                          
                          if (num_mutant > 0 && length(cleaned_column) > num_mutant) {
                            ci <- binom.test(num_mutant, length(cleaned_column), conf.level = 0.95)$conf.int
                            return(c(Prevalence = prevalence * 100,
                                     LowerCI = ci[1] * 100,
                                     UpperCI = ci[2] * 100,
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          } else {
                            return(c(Prevalence = prevalence * 100,
                                     LowerCI = NA,
                                     UpperCI = NA,
                                     NumMutant = num_mutant,
                                     NumWildtype = num_wildtype))
                          }
                        }
                        
                        # Calculate prevalence per Woreda and per marker
                        results_list <- list()
                        
                        woredas <- unique(data1$Woreda)
                        
                        for (woreda in woredas) {
                          woreda_data <- data1[data1$Woreda == woreda, ]
                          
                          for (marker in markers) {
                            marker_column <- woreda_data[[marker]]
                            marker_results <- calculate_prevalence_with_ci(marker_column)
                            
                            results_df <- data.frame(Woreda = woreda,
                                                     Marker = marker,
                                                     Family = "k13",
                                                     Prevalence = marker_results["Prevalence"],
                                                     LowerCI = marker_results["LowerCI"],
                                                     UpperCI = marker_results["UpperCI"],
                                                     NumMutant = marker_results["NumMutant"],
                                                     NumWildtype = marker_results["NumWildtype"])
                            
                            results_list[[paste(woreda, marker, sep = "_")]] <- results_df
                          }
                        }
                        
                        # Combine all results into a single dataframe
                        results_df_final <- do.call(rbind, results_list)
                        
                        # Save to CSV
                        write.csv(results_df_final, "k13_prevalence_by_site.csv", row.names = FALSE)
                        
                        # Optional: Plot example for one Woreda or facet by Woreda
                        ggplot(results_df_final, aes(x = Marker, y = Prevalence, fill = Marker)) +
                          geom_bar(stat = "identity") +
                          geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
                          scale_fill_manual(values = marker_color_map) +
                          facet_wrap(~ Woreda, scales = "free_y") +
                          theme_minimal() +
                          labs(title = "K13 Mutation Prevalence by Site Level",
                               x = "Markers",
                               y = "Prevalence of Mutants (%)",
                               fill = "Marker") +  # <- this labels the legend
                          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                                legend.title = element_text(size = 12),
                                legend.text = element_text(size = 10),
                                legend.position = "right")  # <- shows the color key on the right
                        
                        
                        ################################
                        ################################
                        #metadrcoimergedata <- ("metadrcoimergedata.csv")
                        ################################
                        ################################
                        
                        # COI over all 
                        
                        getwd()
                        
                        metadrcoimergedata <- read.csv("metadrcoimergedata.csv")
                        
                        ggplot(metadrcoimergedata, aes(coi)) +
                          geom_histogram(binwidth = 0.5, fill = "#008080") +
                          labs(y = "Number of Samples", x = "COI") +
                          scale_x_continuous(breaks = seq(floor(min(metadrcoimergedata$coi)), ceiling(max(metadrcoimergedata$coi)), by = 1))+
                          
                          ggtitle("COI distribution based maxhapcount")
                        
                        
                        # Create the histogram plot by site
                        ggplot(metadrcoimergedata, aes(coi)) +
                          geom_histogram(binwidth = 0.5, fill = "#008080", color = "black") +
                          labs(y = "Number of Samples", x = "COI") +
                          scale_x_continuous(
                            breaks = seq(floor(min(metadrcoimergedata$coi)), ceiling(max(metadrcoimergedata$coi)), by = 1)
                          ) +
                          ggtitle("COI Distribution by Site") +
                          facet_wrap(~ site)  # Facet by site
                        
                        
                        
                        #############
                        # Site level COI
                        ############
                        Bokre.Site = c("#e6194b",  "#f58231", "#ffe119",   "#911eb4", 
                                                "#f032e6", "#4363d8","#46f0f0","#008080", "#3cb44b") 
                                                
                        ggplot(metadrcoimergedata, aes(x = site, y = coi, fill= site)) + 
                          geom_violin(trim=FALSE)+
                          #scale_fill_manual(values=site)+
                          scale_fill_manual(values = c("#e6194b",  "#f58231", "#ffe119",   "#911eb4", 
                                                                "#f032e6", "#4363d8","#46f0f0","#008080", "#3cb44b")) +
                                                                  ylab("COI") +
                          theme(axis.text.x=element_text(angle=45,hjust=1))+
                          ggtitle("Complexity of infection Regional Level")
                        
                        
                        
                        # Load necessary library
                        library(ggplot2)
                        
                        # Read your data
                        metadrcoimergedata <- read.csv("metadrcoimergedata.csv")
                        
                        # Create a boxplot with dots showing values
                        ggplot(metadrcoimergedata, aes(x = site, y = coi, fill = site)) + 
                          geom_boxplot(outlier.shape = NA, color = "black", alpha = 0.6) +  # Boxplot without outliers (dots will be used)
                          geom_jitter(width = 0.2, color = "black", size = 1.5) +  # Add dots for individual values
                          scale_fill_manual(values = c("#e6194b", "#f58231", "#ffe119", "#911eb4", 
                                                                "#f032e6", "#4363d8", "#46f0f0", "#008080", "#3cb44b")) +
                                                                  ylab("COI") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                          ggtitle("Complexity of Infection at Regional Level")
                        
                        
                        # Load necessary library
                        library(ggplot2)
                        
                        # Read your data
                        metadrcoimergedata <- read.csv("metadrcoimergedata.csv")
                        
                        # Create a boxplot without dots
                        ggplot(metadrcoimergedata, aes(x = site, y = coi, fill = site)) + 
                          geom_boxplot(color = "black", alpha = 0.6) +  # Boxplot without individual dots
                          scale_fill_manual(values = c("#e6194b", "#f58231", "#ffe119", "#911eb4", 
                                                                "#f032e6", "#4363d8", "#46f0f0", "#008080", "#3cb44b")) +
                                                                  ylab("COI") +
                          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                          ggtitle("Complexity of Infection at Regional Level")
                        
                        # Read the data
                        
                        
                        #####################
                        #####################
                        #COI
                        #####################
                        #####################
                        
                        library(McCOILR)
                        library(tidyverse)
                        
                        #example dataset
                        
                        data1= read.table("COIFROMATED05.txt", sep="", head=T)
                        
                        # Modify values between 0 and 1 to 0.5
                        data1[] <- lapply(data1, function(x) {
                          # Replace values between 0 and 1 with 0.5
                          x[x > 0 & x < 1] <- 0.5
                          return(x)
                        })
                        
                        # Save the modified data to a new text file
                        write.table(data1, "coiformated1.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
                        
                        # Print a message to confirm the file has been saved
                        cat("The modified data has been saved as 'coiformated1.txt'.")
                        
                        
                        # Overall COI
                        
                        data1= read.table("coiformated1.txt", sep="", head=T)
                        data2=data1[,-1]
                        rownames(data2)=data1[,1]
                        set.seed(2024)
                        
                        out_cat <- McCOIL_categorical(data2,
                                                      maxCOI = 25, 
                                                      threshold_ind = 20, 
                                                      threshold_site = 20,
                                                      totalrun = 1000, 
                                                      burnin = 100, 
                                                      M0 = 15, 
                                                      e1 = 0.05, 
                                                      e2 = 0.05,
                                                      err_method = 3,
                                                      path = "COI", 
                                                      output="Overall_COI.tsv")
                        
                        
                        data2 <- read.csv("Overall_COI.csv")
                        # Have a look at the categorical output COI distribution
                        hist(
                          as.numeric(
                            as.character(
                              data2$median
                            )),
                          main = "Overall COI",  
                          xlab="Median COI", col= "darkred")
                        
                        hist(
                          as.numeric(
                            as.character(
                              data2$mean
                            )),
                          main = "Overall Mean COI",  
                          xlab="Mean COI", col="#008080")
                        
                        
                        ggplot(data2, aes(mean)) +
                          geom_histogram(binwidth = 0.5, fill = "#008080") +
                          labs(y = "Number of Samples", x = "COI") +
                          scale_x_continuous(breaks = seq(floor(min( data2$mean)), ceiling(max( data2$mean)), by = 1))+
                          
                          ggtitle("COI distribution based THEREALMcCOIL")
                        
                        
                        ##
                        # More analysis coming 
                        ##
                        
                        
                        