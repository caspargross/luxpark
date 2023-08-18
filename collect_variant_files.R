# Imoprt files
library(data.table)
library(tidyverse)


dtl_pb <- list()
# Import PacBio Deepvariant
samples <- str_replace(list.files("../PacBio", pattern="^Sample_*"), "Sample_", "")
for (sample in samples[!str_detect(samples, "_2$")]) {
  dtl_pb[[sample]] <- read_tsv(paste0("../PacBio/Sample_", sample,"/", sample, ".dv.tsv"),
                            col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "GQ", "DP", "AD", "VAF", "PL"),
                            col_types = c("cdccdccddccc")) %>%
    mutate(sample = str_replace(sample, "_2", "")) %>% 
    mutate(sequencing_method = "SequelII") %>%
    mutate(variant_caller = "Deepvariant PB") %>%
    mutate(batch = "B1") %>%
    mutate(replicate = "1") 
}
for (sample in samples[str_detect(samples, "_2$")]) {
  dtl_pb[[sample]] <- read_tsv(paste0("../PacBio/Sample_", sample,"/", sample, ".dv.tsv"),
                            col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "GQ", "DP", "AD", "VAF", "PL"),
                            col_types = c("cdccdccddccc")) %>%
    mutate(sample = str_replace(sample, "_2", "")) %>% 
    mutate(sequencing_method = "SequelII") %>%
    mutate(variant_caller = "Deepvariant PB") %>%
    mutate(batch = "B1") %>%
    mutate(replicate = "2")
}

dt_pb <- bind_rows(dtl_pb)
samples_pacbio <- read_csv("sample_names_pacbio.csv")
dt_pb <- left_join(dt_pb, samples_pacbio %>% select(sample = "Name intern", sample_extern = "Name extern"))

# Import Nanopore Medaka
files <- Sys.glob("../Sample_*/demux_analysis/*/round_1.tsv")
dtl_ont <- list()

for (file in files) {
      
    sample <- str_match(file, "20010Ea00[1,2]_(.*)/round")[2]
    print(sample)
    dtl_ont[[file]] <- read_tsv(file,
             col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "GQ"),
             col_types = c("cdccdccd")) %>%
    mutate(sequencing_method = ifelse(str_detect(file, "20010Ea002"),"MinIon LowCov", "MinIon")) %>%
    mutate(variant_caller = "Medaka") %>%
    mutate(batch = "B1") %>%
    mutate(replicate = ifelse(str_detect(sample, "_2"),"2", "1")) %>%
    mutate(sample = str_replace(sample, "_2", "")) 
}

dt_ont <- bind_rows(dtl_ont)
samples_minion <- read_csv("sample_names_minion.csv")
dt_ont <- left_join(dt_ont, samples_minion %>% select(sample = "Name intern", sample_extern = "Name extern"))


# Import Nanopore Pepper Deepvariant
files <- Sys.glob("../variants_pepper/*/pepper_variants.tsv")
dtl_ontdv <- list()

for (file2 in files) {
  sample <- str_match(file2, "20010Ea00[1,2]_(.*)/pepper_variants.tsv")[2]
  print(paste(file2, str_detect(file2, "20010Ea002"))) 
  print(sample)
  method <- ifelse(str_detect(file2, "20010Ea002"), "MinIon LowCov", "MinIon")
  print(method)
  dtl_ontdv[[file2]] <- read_tsv(file2,
                                col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "GQ", "DP", "AD", "VAF", "PL"),
                                col_types = c("cdccdccddccc")) %>%
    mutate(variant_caller = "DeepVariant PEPPER") %>%
    mutate(batch = "B1") %>%
    mutate(sequencing_method = method) %>%
    mutate(replicate = ifelse(str_detect(sample, "_2"),"2", "1")) %>%
    mutate(sample = str_replace(sample, "_2", "")) 
}

dt_ontdv <- bind_rows(dtl_ontdv)
samples_minion <- read_csv("sample_names_minion.csv")
dt_ontdv <- left_join(dt_ontdv, samples_minion %>% select(sample = "Name intern", sample_extern = "Name extern"))


## TODO
dt <- bind_rows(dt_ont, dt_pb %>% filter(FILTER == "PASS"), dt_ontdv %>% filter(FILTER == "PASS"))


write_tsv(dt, "variants_all_GBA_PASS.tsv")

# plot  
p <- ggplot(dt, aes(x=POS, col = ALT, fill = ALT)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(limits = c(155232501,155241415)) +
  theme_classic()
p

ggsave("GBA.png", width=20)


dtl_pbsv <- list()
# Import PacBio PBSV
# %CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%AD\t%DP\t%SAC\t%CN
samples <- str_replace(list.files("../PacBio_B2/demux_m64074_201210_132215/", pattern="^Sample_*"), "Sample_", "")
samples <- samples[!samples[] == "unknown"]
for (sample in samples[!str_detect(samples, "_2$")]) {
  print(sample)
  dtl_pbsv[[sample]] <- read_tsv(paste0("../PacBio_B2/demux_m64074_201210_132215/Sample_", sample,"/", sample, ".pbsv.tsv"),
                            col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "AD", "DP", "SAC", "CN"),
                            col_types = c("cdccdccddcc")) %>%
    mutate(sample = str_replace(sample, "_2", "")) %>% 
    mutate(sequencing_method = "SequelII") %>%
    mutate(variant_caller = "Deepvariant PB") %>%
    mutate(batch = "B1") %>%
    mutate(replicate = "1") 
}
for (sample in samples[str_detect(samples, "_2$")]) {
  print(sample)
  dtl_pbsv[[sample]] <- read_tsv(paste0("../PacBio_B2/demux_m64074_201210_132215/Sample_", sample,"/", sample, ".pbsv.tsv"),
                            col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "AD", "DP", "SAC", "CN"),
                            col_types = c("cdccdccddcc")) %>%
    mutate(sample = str_replace(sample, "_2", "")) %>% 
    mutate(sequencing_method = "SequelII") %>%
    mutate(variant_caller = "Deepvariant PB") %>%
    mutate(batch = "B1") %>%
    mutate(replicate = "2")
}
samples <- str_replace(list.files("../PacBio_B3/demux_20010E_Batch3/", pattern="^Sample_*"), "Sample_", "")
samples <- samples[!samples[] == "unknown"]
for (sample in samples[!str_detect(samples, "_2$")]) {
  print(sample)
  dtl_pbsv[[sample]] <- read_tsv(paste0("../PacBio_B3/demux_20010E_Batch3/Sample_", sample,"/", sample, ".pbsv.tsv"),
                            col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "AD", "DP", "SAC", "CN"),
                            col_types = c("cdccdccddcc")) %>%
    mutate(sample = str_replace(sample, "_2", "")) %>% 
    mutate(sequencing_method = "SequelII") %>%
    mutate(variant_caller = "Deepvariant PB") %>%
    mutate(batch = "B2") %>%
    mutate(replicate = "1") 
}
for (sample in samples[str_detect(samples, "_2$")]) {
  print(sample)
  dtl_pbsv[[sample]] <- read_tsv(paste0("../PacBio_B3/demux_20010E_Batch3/Sample_", sample,"/", sample, ".pbsv.tsv"),
                            col_names = c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "GT", "AD", "DP", "SAC", "CN"),
                            col_types = c("cdccdccddcc")) %>%
    mutate(sample = str_replace(sample, "_2", "")) %>% 
    mutate(sequencing_method = "SequelII") %>%
    mutate(variant_caller = "Deepvariant PB") %>%
    mutate(batch = "B2") %>%
    mutate(replicate = "2")
}

dt_pbsv <- bind_rows(dtl_pbsv)


samples_pacbio <- read_csv("../analysis/sample_names_pacbio.csv")
dt_pbsv <- left_join(dt_pbsv, samples_pacbio %>% select(sample = "Name intern", sample_extern = "Name extern"))

fwrite(dt_pbsv, "summary_table_sv.tsv", sep="\t")
