library(ggplot2)
library(dplyr)

theme_set(theme_bw())

args <- commandArgs(trailingOnly = TRUE)

tsv <- args[1]
workdir <- args[2]
prefix <- args[3]
datatype <- args[4]

setwd(workdir)

df <- read.table(tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

totalreads_int <- as.integer(sum(df$num_reads))

# adding commas for every thousanths place present in totalreads var
totalreads <- format(totalreads_int, big.mark = ",", scientific = FALSE)

# trim any leading spaces from 3rd column in dataframe
df$percent_reads <- sub("^\\s+", "", df$percent_reads)

# convert 3rd column to double
df$percent_reads <- as.double(df$percent_reads)

# grab the top 20 rows by total_species_count
unique_species <- length(unique(df$species))
if (unique_species > 20) {

  top_species <- df %>%
  slice_max(n = 15, order_by = num_reads)

  # if we want to grab the top 20 species we want to make sure we also grab all ambiguity entries of each species, even if the read count for one ambiguity level is much lower than the other
  for (i in 1:nrow(df)) {

    row <- df[i, ]

    if (row$species %in% top_species$species) {

      if (exists("top_hits")) {

        top_hits <- rbind(top_hits, row)

      } else {

        top_hits <- row

      }

    }
  
  }

} else {

  top_hits <- df

}

top_hits_mod <- top_hits %>%
  arrange(species)

### This block will calculate the number and percent of "other" species classified by the "minimap2" algorithm and add them as a new row
total_percent <- sum(top_hits_mod$percent_reads)

if (total_percent < 99.99) {

  # reference original df to grab read count of species (not just top 20)
  read_count_all <- sum(df$num_reads)

  percent_other <- 100 - total_percent
  count_other <- read_count_all - sum(top_hits_mod$num_reads)

  new_row_other <- data.frame(
  species = "other",
    num_reads = count_other,
    percent_reads = percent_other,
    ambiguity = "none"
  )
  top_hits_mod <- rbind(top_hits_mod, new_row_other)

}

color_list <- c(
  "#FF7666" # Lighter Red
)

top_hits_mod$color <- NA

#setting up a separate dataframe that will make it easier to plot with ggplot
ggplot_df <- top_hits_mod

ggplot_df$color[ggplot_df$ambiguity == "none"] <- "#FF7666"
ggplot_df$color[ggplot_df$ambiguity == "single_genome"] <- "#ADD8E6"
ggplot_df$color[ggplot_df$ambiguity == "multi_genome"] <- "#D8BFD8"

ggplot_df$total_species_count <- NA
ggplot_df$total_species_percent <- NA
ggplot_df$offset <- NA

# creating samll df to see which species contains the highest read percentage regardless of ambiguity
total_count_by_species <- aggregate(num_reads ~ species, data = top_hits_mod, sum)
total_percent_by_species <- aggregate(percent_reads ~ species, data = top_hits_mod, sum)

# grabbing highest percent value to adjust data label positions in plot. The highest percent value will ultimately determine the scale of the y-axis, so our data label positions will be based on that
max_percent = max(total_percent_by_species$percent_reads, na.rm = TRUE)

# size factor is used to estimate the height of a single data label - the number 15 was determined via trial and error
size_factor = max_percent / 15

# this for loop will determine the y-axis position of each data label to ensure there is no overlap
for (i in 1:nrow(ggplot_df)) {

  row <- ggplot_df[i, ]

  ggplot_df[i, "total_species_count"] <- total_count_by_species$num_reads[total_count_by_species$species == row$species]

  ggplot_df[i, "total_species_percent"] <- total_percent_by_species$percent_reads[total_percent_by_species$species == row$species]

  if (row$ambiguity == "none") {

    ggplot_df[i, "offset"] <- (total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - ( total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - row$percent_reads)) / 2

  } else if (row$ambiguity == "single_genome") {

    none_ambiguity_percent = ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "none"]

    vertical_pos = (row$percent_reads / 2) + none_ambiguity_percent

    if (length(none_ambiguity_percent) == 0) {

      vertical_pos = vertical_pos

    } else if (none_ambiguity_percent >= size_factor * 2) {

      vertical_pos = vertical_pos

    } else if (none_ambiguity_percent < size_factor | row$percent_reads < size_factor) {

      vertical_pos = vertical_pos + size_factor

    }

    ggplot_df[i, "offset"] <- vertical_pos

  } else {

    none_ambiguity_percent = ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "none"]
    single_ambiguity_percent = ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "single_genome"]
    single_ambiguity_offset = ggplot_df$offset[ggplot_df$species == row$species & ggplot_df$ambiguity == "single_genome"]

    vertical_pos = total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - (row$percent_reads / 2)

    # determines if 'single_genome' and 'none' ambiguity classifiers exists for this species; this is required to prevent a 'length of zero' error in the subsequent if statements
    if (length(single_ambiguity_offset) == 0 & length(none_ambiguity_percent) == 0) {

      vertical_pos = vertical_pos

    } else if (length(single_ambiguity_offset) == 0 & length(none_ambiguity_percent) != 0) {

      if (none_ambiguity_percent >= size_factor) {

        vertical_pos = vertical_pos
      
      } else {

        vertical_pos = none_ambiguity_percent + size_factor

      }

    } else if (length(single_ambiguity_offset) != 0 & length(none_ambiguity_percent) == 0) {

        if (single_ambiguity_offset + size_factor <= vertical_pos) {

          vertical_pos = vertical_pos
        
        } else {

          vertical_pos = single_ambiguity_offset + size_factor

        }

    } else if (single_ambiguity_offset + size_factor <= vertical_pos) {

      vertical_pos = vertical_pos
    
    } else {

      vertical_pos = single_ambiguity_offset + size_factor

    }

    ggplot_df[i, "offset"] <- vertical_pos

  }

}

# force ggplot2 to use the order of species as they appear in ggplot_df by converting species column to a factor with levels set in the desired order
ggplot_df$species <- factor(ggplot_df$species, levels = unique(ggplot_df$species))
ggplot_df$ambiguity <- factor(ggplot_df$ambiguity, levels = rev(unique(ggplot_df$ambiguity)))


plot_1 <- ggplot(ggplot_df, aes(x = species, group = ambiguity, y = percent_reads)) +
  geom_col(aes(fill = ambiguity), color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  scale_fill_manual(values = setNames(ggplot_df$color, ggplot_df$ambiguity)) +
  labs(
    title = paste(prefix, datatype, "Taxonomic Distribution"),
    subtitle = paste0("Total ", datatype, "s: ", totalreads),
    x = "",
    y = paste0("Percent Mapped")
  ) +
  geom_label(
    data = ggplot_df,
    size = 2.8,
    aes(label=paste0(sprintf("%.2f", percent_reads), "%", " | ", num_reads), y = offset, fill = factor(ambiguity)),
    label.padding = unit(0.15, "lines"),
    show.legend  = FALSE
  )

ggsave("plot.png", width = 16, height = 6, dpi = 150)
