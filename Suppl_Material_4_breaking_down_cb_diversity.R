#### TÍTULO: Procesamiento de amplicones 16S, restauración de manglares en Celestún  ###
#### Librerías requeridas: qiime2R, phyloseq, tidyverse, ggplot2, dplyr, vegan ###
#### Autores: Alejandro Ávila, Daniel Vázquez, Daniel Esguerra ####
#### Creado el: ? ####
#### Última modifcación: 22 de marzo de 2025 ####

#Vamos a tomar como referencia el codigo de Daniel Esguerra para hacer el analisis de los
#resultados de quiime para manglares en estatus de restauración de Celestun 

#### Librerias ####
library(qiime2R)
library(phyloseq)
library(picante)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(vegan)
library(ade4)
library(microbiome)
library(plotly)
library(FSA)
library(ggpubr)
library(VIM)
library(MetBrewer)
library(rstatix)

# install.packages("processx")  # Required dependency
# devtools::install_github("plotly/orca")

#### Preparar datos ####

#establecemos el directorio de trabajo
#setwd("~/maestria/Proyecto Manglares Celestun/Procesamiento/procesamiento-211")

### Abrimos datos como un phyloseq object ###
physeq<-qza_to_phyloseq(
  features="./table.qza",
  tree="./rooted-tree.qza",
  taxonomy="./taxonomy5.qza",
  metadata = "./mapa_restauracionC.txt"
)

### Filtros de datos... ####
#Establecemos un seed para que los datos sean reproducibles
set.seed(123)

#remove those low-quality samples by means of the sampling depth obtained by
#DADA2 filter
physeq_filt <- rarefy_even_depth(physeq, sample.size = 10117) #Alex: 16828 #>1000 reads: 10117

# 23 samples removedbecause they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 
#   
#   R1.1.1R1.2.0R2.1.1R2.1.2R2.3.2
# ...
# 3971 OTUs were removed because they are no longer (3936 OTUs with 17111)
# present in any sample after random subsampling

# 2) Data transformation and cleaning

#we remove mock samples
# defined mixtures of known microbial cells, DNA, or RNA, either physically 
# created or computationally simulated, used as ground truth to validate and 
# optimize analysis methods by providing a known composition to compare against 
# experimental results
physeq_filt <- prune_samples(sample_names(physeq_filt) != "Mock", physeq_filt) # Esto debe filtrar muestras que no queremos

#Vamos a filtrar elementos no asignados a la taxonomia en Kingdom 
physeq_filt<-physeq_filt %>% subset_taxa(Kingdom != "Unassigned")

#Filtar las Eucariontes y mitocondrias, lo cual podría estar relacionado con contaminación
#de DNA y problemas con las bases de datos
filtphyla<-c("Mitochondria", "Eukaryota", "Chloroplast")
physeq_filt <- subset_taxa(physeq_filt , !Kingdom %in% filtphyla)
physeq_filt <- subset_taxa(physeq_filt , !Phylum %in% filtphyla)
physeq_filt <- subset_taxa(physeq_filt , !Class %in% filtphyla)
physeq_filt <- subset_taxa(physeq_filt , !Order %in% filtphyla)
physeq_filt <- subset_taxa(physeq_filt , !Family %in% filtphyla)

#filtrar muestras con un %menor a 10,000 nucletidos despues de DADA2
#y que tienen una composición taxonómica rara en el barplot de qiime2

#creamos el vector con las muestras a remover
ss<-c("R2.1.1", "R1.1.1", "R4.3.1",
      "R4.2.3", "R4.1.3", "R5.1.1",
      "R6.1.1", "R3.2.1", "R4.2.1",
      "R5.3.1", "R1.2.0", "R5.2.3",
      "R0.3.3", "R2.1.2", "R7.3.2",
      "R7.4.1", "R7.4.2", "R7.4.3") #the last three correspond to z4

physeq_filt <- prune_samples(!(sample_names(physeq_filt) %in% ss), physeq_filt)

#Como queda antes y despues de filtros 
tax_antes<-physeq@tax_table@.Data
tax_desp<-physeq_filt@tax_table@.Data

diferencia<-length(tax_antes)-length(tax_desp)
diferencia
#Se fueron 28252 ASV que no se asignaron taxonomicamente (28119 a 17111)

#rm(tax_antes,tax_desp,diferencia)

#we extract the taxonomy table to add the assignations per placement
tax_df <- as.data.frame(tax_table(physeq_filt))
write.csv(tax_df, "./tax_table.csv")

#import the modified table
tax_df_cb <- (read.csv("./tax_table_mod.csv"))
asvs <- tax_df_cb[,1]
tax_df_cb <- tax_df_cb[,c(-1)]
rownames(tax_df_cb) <- asvs
#assign it to physeq obj
tax_table(physeq_filt) <- as.matrix(tax_df_cb)

#Vamos a dejar los metadatos a la mano para facilitar el manejo de datos 
metadata_ <- as.data.frame(physeq_filt@sam_data)
sample_id <- rownames(metadata_)
metadata_ <- cbind(sample_id, metadata_)
rownames(metadata_) <- NULL
rm(sample_id)

#traducimos las etiquetas al inglés para mayor consistencia
metadata_$Description[metadata_$Description == "Conservado"] <- "Conserved"
metadata_$Description[metadata_$Description == "Restaurado"] <- "Restored"
metadata_$Description[metadata_$Description == "Degradado"] <- "Degraded"

#y añadimos los rangos en cm a los Z layers
metadata_$Z_Depth[metadata_$Z_Depth == "Z1"] <- "Z1 (00-15 cm)"
metadata_$Z_Depth[metadata_$Z_Depth == "Z2"] <- "Z2 (15-30 cm)"
metadata_$Z_Depth[metadata_$Z_Depth == "Z3"] <- "Z3 (30-50 cm)"
metadata_$Z_Depth[metadata_$Z_Depth == "Z4"] <- "Z4 (50-100 cm)"

#creamos un gráfico de barras apilado de composición con filos >1%
physeq_dec <- transform_sample_counts(physeq_filt, function(x){x / sum(x)})
physeq_porcen <- transform_sample_counts(physeq_dec, function(x){x * 100})

#
otu_counts <- physeq_filt@otu_table@.Data
otu_porcen <- physeq_dec@otu_table@.Data

#extract cb counts
cbs_counts <- otu_counts[rownames(otu_desp) %in% c('53d9c860ca705176b1f2861dac2ddecc', #yucatan.7
                               '8e8a54e4274b6f0578ffc60725b2d149', #yucatan.8
                               'ce101702ae7b495905154e400eca0f8f', #yucatan.10
                               '70d706b29038d62a218820aedca7a388'), , drop = FALSE] #yucatan.13

cbs_porcen <- otu_porcen[rownames(otu_porcen) %in% c('53d9c860ca705176b1f2861dac2ddecc', #yucatan.7
                        '8e8a54e4274b6f0578ffc60725b2d149', #yucatan.8
                        'ce101702ae7b495905154e400eca0f8f', #yucatan.10
                        '70d706b29038d62a218820aedca7a388'), , drop = FALSE] #yucatan.13

read_counts <- sample_sums(physeq_filt)

#merging taxonomy and metadata

#taxonomy
cbs_asvs <- rownames(cbs_counts)
#subset tax table
cbs_tax <- tax_df_cb[rownames(tax_df_cb) %in% cbs_asvs, , drop = FALSE]
#merge with taxonomy
#tax + counts
cbs_counts <- cbind(as.data.frame(cbs_counts), cbs_asvs, cbs_tax)
rownames(cbs_counts) <- NULL

#tax + per
cbs_porcen <- cbind(as.data.frame(cbs_porcen), cbs_asvs, cbs_tax)
rownames(cbs_porcen) <- NULL

#make matrices long
#counts
cb_counts_long <- pivot_longer(cbs_counts, cols = c(1:60), names_to = 'sample_id', values_to = 'count')
#perc
cbs_porcen_long <- pivot_longer(cbs_porcen, cols = c(1:60), names_to = 'sample_id', values_to = 'percentage')

#and join metadata
cb_counts_full <- inner_join(cb_counts_long, metadata_,by = 'sample_id')
cb_porcen_full <- inner_join(cbs_porcen_long, metadata_, by = 'sample_id')

desc <- c("Conserved", "Degraded", "Restored")
cb_counts_desc <- aggregate(cb_counts_full, count ~ Description + Genus, FUN = sum)
cb_percen_desc <- cb_counts_desc
cb_percen_desc$count[1] <- round((cb_counts_desc$count[1] / read_counts[c('R5.2.1')]) *100, digits = 4)
cb_percen_desc$count[3] <- round((cb_counts_desc$count[3]) / sum(read_counts[c('R4.3.2', 'R4.1.2', 'R3.3.3', 'R4.1.2d', 'R3.3.1')]) *100, digits = 4)
cb_percen_desc$count[6] <- round((cb_counts_desc$count[6]) / sum(read_counts[c('R4.3.2', 'R4.1.2', 'R3.3.3', 'R4.1.2d', 'R3.3.1')]) *100, digits = 4)

#Modify names
cb_percen_desc$Genus[1:3] <- "Cluster I"
cb_percen_desc$Genus[4:6] <- "Cluster VI"

cb_prevalence <- ggplot(cb_percen_desc, aes(y = Description, x = count, fill = Genus)) +
  geom_bar(position = "dodge", stat= "identity") +
  scale_fill_manual(values = c("#4ed322","#cf9f0e")) +
  labs(#title = "Phyla-level abundance composition", #título general
    x = "Relative abundance", #título para eje "y"
    y = "Status",
    fill = "Genus-level clades") + #usamos título en eje "x" para la posición geográfica
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(vjust = 0.5, hjust=1)) +
  guides(fill = guide_legend(ncol = 2)) #ordenamos etiquetas de valores (phyla) en 5 columnas

cb_prevalence
ggsave("./cables_diversity/cb_prevalence_celestun.png", cb_prevalence, 
       width = 27, height = 20, units = "cm", dpi = 600)

#let's make an upset plot for the diversity per location (worldwide)
library(UpSetR)
library(ComplexUpset)

richness_location <- read.csv("matrix_locations.csv")
richness_per_site <- upset(
  richness_location,
  c("Cluster.V",
     "Cluster.IV",
     "Cluster.VI",
     "Cluster.III",
     "Cluster.I",
     "Cluster.II",
     "Unclassified_cb"),
  queries = list(
    upset_query(set='Cluster.V', fill='#029a9a'),
    upset_query(set='Cluster.IV', fill='#a115a1'),
    upset_query(set='Cluster.VI', fill='#cf9f0e'),
    upset_query(set='Cluster.III', fill='#75751a'),
    upset_query(set='Cluster.I', fill='#4ed322'),
    upset_query(set='Cluster.II', fill='#099d6b'),
    upset_query(set='Unclassified_cb', fill='#666666')
  ),
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        bar_number_threshold=1,  # show all numbers on top of bars
        width=0.5,   # reduce width of the bars
      )
      # add some space on the top of the bars
      + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
      + theme(
        # hide grid lines
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # show axis lines
        axis.line=element_line(colour='black')
      )
    )
  ),
  stripes=upset_stripes(
    geom=geom_segment(linewidth=12),  # make the stripes larger
    colors=c('grey95', 'white')
  ),
  # to prevent connectors from getting the colorured
  # use `fill` instead of `color`, together with `shape='circle filled'`
  matrix=intersection_matrix(
    geom=geom_point(
      shape='circle filled',
      size=3.5,
      stroke=0.45
    )
  ),
  set_sizes=(
    upset_set_size(geom=geom_bar(width=0.4))
    + theme(
      axis.line.x=element_line(colour='black'),
      axis.ticks.x=element_line()
    )
  ),
      sort_sets='descending',
      sort_intersections = 'descending'
      )

richness_per_site
ggsave("./cables_diversity/locations_breakdown.png", richness_per_site, 
       width = 38, height = 24, units = "cm", dpi = 600)

#bubble plot
library(dplyr)

locations_habitats <- read.csv("./locations_habitats.csv")
habitats_bubble <- locations_habitats %>%
  group_by(main_habitat, assigned) %>%
  summarise(
    locations = n_distinct(location_num),
    .groups = "drop"
  )


habitats_location_bubble <- ggplot(habitats_bubble,
                                   aes(x=assigned, y=main_habitat, 
                                       size = locations, color=assigned)) +
  geom_point(alpha = 0.7) +
  scale_size(range = c(1, 15), name=("Num. of Locations")) +
  scale_color_manual(values = c('#4ed322', '#099d6b', '#75751a',
                                '#a115a1', '#029a9a', '#cf9f0e',
                                '#666666')) +
  labs(
    y = "Main habitat", #título para eje "y"
    x = "Cable bacteria clusters") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(vjust = 1, hjust=1, angle = 45)) +
  guides(color = "none") #ordenamos etiquetas de valores (phyla) en 5 columnas

habitats_location_bubble
ggsave("./cables_diversity/habitats_per_location.png", 
       habitats_location_bubble, width = 27, height = 24, 
       units = "cm", dpi = 600)

#richness

#per location
richness_per_location <- locations_habitats %>%
  group_by(main_habitat, location_num) %>%
  summarise(
    richness = n_distinct(assigned),
    .groups = "drop"
  )

#per habitat
richness_per_habitat <- richness_per_location %>%
  group_by(main_habitat, richness) %>%
  summarise(
    amount_locations = n_distinct(location_num),
    .groups = "drop"
  )

richness_habitat_bubbl <- ggplot(richness_per_habitat,
                                 aes(x = amount_locations, y = main_habitat, 
                                     size = amount_locations,
                                     color = richness)) +
  geom_point(alpha = 0.85) +
  scale_size(range = c(1, 20), name=("Number of Locations")) +
  scale_color_met_c("Johnson",direction = -1) +
  labs(
    y = "Main Habitat", #título para eje "y"
    x = "Num. Locations",
    color = "Richness") +
  theme_bw() +
  theme(legend.position = "bottom") 

richness_habitat_bubbl
ggsave("./cables_diversity/richness_per_habitat.png", richness_habitat_bubbl, width = 27, height = 24, units = "cm", dpi = 600)
