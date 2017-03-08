# Libraries ---------------------------------------------------------------
library(dplyr)
library(magrittr)
library(stringr)
library(reshape2)
library(data.table)
load('Human.info.Rdata')
load('Mouse.info.Rdata')
# Generate (Human)----------------------------------------------------------------
hum.path <-
  read.csv(
    'Human_pathlist.txt',
    header = F,
    sep = "\t",
    stringsAsFactors = F
  )

hum.path <-
  hum.path %>% 
  mutate(index = seq(1, nrow(hum.path))) %>% 
  as.data.table() #hum.path, add a index for merge
hum.path.name <-
  hum.path[!str_detect(V2, '\\w')] %>%
  mutate(V2 = NULL, V3 = NULL) # Separate path name
hum.path.ens <-
  hum.path[V2 == 'Ensembl'] %>%
  left_join(hum.info, by = c('V3' = 'Ensembl.hum')) %>%
  mutate(Symbol.hum = coalesce(Symbol.hum, V1),
         Ensembl.hum = V3,
         V1 = as.character(NA),
         V3 = NULL,
         V2 = NULL) #path with ensembl database
hum.path.ent <-
  hum.path[V2 == 'EntrezGene'] %>%
  mutate(V3 = as.integer(V3)) %>%
  left_join(hum.info, by = c('V3' = 'Entrez.hum')) %>%
  mutate(Symbol.hum = coalesce(Symbol.hum, V1),
         Entrez.hum = V3,
         V1 = as.character(NA),
         V2 = NULL,
         V3 = NULL) # path with entrez database
hum.path.uni <-
  hum.path[V2 == 'Uniprot'] %>%
  left_join(hum.info, by = c('V3' = 'Accession.hum')) %>%
  mutate(
    Symbol.hum = coalesce(Symbol.hum, V1),
    Accession.hum = V3,
    V1 = as.character(NA),
    V2 = NULL,
    V3 = NULL) # path with uniprot database

hum.path <- full_join(hum.path.ens,hum.path.ent) %>% 
  full_join(hum.path.uni) %>% 
  full_join(hum.path.name) %>% 
  arrange(index) %>% 
  select(pathname.hum = V1,
         Ensembl.hum,
         Entrez.hum,
         Symbol.hum,
         Accession.hum)

hum.path[nrow(hum.path) + 1, 1] <- 'Name=stop'
hum.pathname.index <-
  which(str_detect(hum.path[, 1], 'Name=')) #Locate the pathway name.
hum.path[, 1] <- hum.path[, 1] %>% str_replace_all('Name=', '')

hum.pathname <- as.character(NULL)
nullpath <- as.integer(NULL)
hum.pathlist <- list(NULL)
for (i in 1:(length(hum.pathname.index) - 1)) {
  if (hum.pathname.index[i] + 1 == hum.pathname.index[i + 1]) {
    cat(hum.path[hum.pathname.index[i], 'pathname.hum'], 'have no proteins. Omitted.\n')
    nullpath <- c(nullpath, i)
    hum.pathlist[[i]] <- NA
    hum.pathname[i] <- NA
    next()
  }
  
  temp <- hum.path[seq(from = hum.pathname.index[i] + 1, to = hum.pathname.index[i + 1] - 1), ] 
  hum.pathname[i] <- hum.path[hum.pathname.index[i], 1]
  temp$pathname.hum <- hum.pathname[i]
  hum.pathlist[[i]] <- temp
} #Separate the hum.path.

names(hum.pathlist) <- hum.pathname 
hum.pathlist[nullpath] <- NULL #Remove NULL hum.path.


save(hum.pathlist, file = 'Human.pathwaylist.Rdata')
# Generate (Mouse)----------------------------------------------------------------
mus.path <-
  read.csv(
    'Mouse_pathlist.txt',
    header = F,
    sep = "\t",
    stringsAsFactors = F
  )

mus.path <-
  mus.path %>% 
  mutate(index = seq(1, nrow(mus.path))) %>% 
  as.data.table() #mus.path, add a index for merge
mus.path.name <-
  mus.path[!str_detect(V2, '\\w')] %>%
  mutate(V2 = NULL, V3 = NULL) # Separate path name
mus.path.ens <-
  mus.path[V2 == 'Ensembl'] %>%
  left_join(mus.info, by = c('V3' = 'Ensembl.mus')) %>%
  mutate(Symbol.mus = coalesce(Symbol.mus, V1),
         Ensembl.mus = V3,
         V1 = as.character(NA),
         V3 = NULL,
         V2 = NULL) #path with ensembl database
mus.path.ent <-
  mus.path[V2 == 'EntrezGene'] %>%
  mutate(V3 = as.integer(V3)) %>%
  left_join(mus.info, by = c('V3' = 'Entrez.mus')) %>%
  mutate(Symbol.mus = coalesce(Symbol.mus, V1),
         Entrez.mus = V3,
         V1 = as.character(NA),
         V2 = NULL,
         V3 = NULL) # path with entrez database
mus.path.uni <-
  mus.path[V2 == 'Uniprot'] %>%
  left_join(mus.info, by = c('V3' = 'Accession.mus')) %>%
  mutate(
    Symbol.mus = coalesce(Symbol.mus, V1),
    Accession.mus = V3,
    V1 = as.character(NA),
    V2 = NULL,
    V3 = NULL) # path with uniprot database

mus.path <- full_join(mus.path.ens,mus.path.ent) %>% 
  full_join(mus.path.uni) %>% 
  full_join(mus.path.name) %>% 
  arrange(index) %>% 
  select(pathname.mus = V1,
         Ensembl.mus,
         Entrez.mus,
         Symbol.mus,
         Accession.mus)

mus.path[nrow(mus.path) + 1, 1] <- 'Name=stop'
mus.pathname.index <-
  which(str_detect(mus.path[, 1], 'Name=')) #Locate the pathway name.
mus.path[, 1] <- mus.path[, 1] %>% str_replace_all('Name=', '')

mus.pathname <- as.character(NULL)
nullpath <- as.integer(NULL)
mus.pathlist <- list(NULL)
for (i in 1:(length(mus.pathname.index) - 1)) {
  if (mus.pathname.index[i] + 1 == mus.pathname.index[i + 1]) {
    cat(mus.path[mus.pathname.index[i], 'pathname.mus'], 'have no proteins. Omitted.\n')
    nullpath <- c(nullpath, i)
    mus.pathlist[[i]] <- NA
    mus.pathname[i] <- NA
    next()
  }
  
  temp <- mus.path[seq(from = mus.pathname.index[i] + 1, to = mus.pathname.index[i + 1] - 1), ] 
  mus.pathname[i] <- mus.path[mus.pathname.index[i], 1]
  temp$pathname.mus <- mus.pathname[i]
  mus.pathlist[[i]] <- temp
} #Separate the mus.path.

names(mus.pathlist) <- mus.pathname 
mus.pathlist[nullpath] <- NULL #Remove NULL mus.path.
save(mus.pathlist, file = 'Mouse.pathwaylist.Rdata')

# allpath information from CPDB ----------------------------------------------------------
str2list <- function(chr){
  as.character(str_split(chr,',',simplify = T))
}

hum.allpath <- fread('CPDB_pathways_genes.tab')
hum.allpath.name <- paste(hum.allpath$pathway,hum.allpath$external_id,hum.allpath$source,sep = '_')
hum.allpath.list <- lapply(hum.allpath$entrez_gene_ids,str2list)
names(hum.allpath.list) <- hum.allpath.name
hum.path.genes.ent <- unlist(hum.allpath.list) %>% unique()
hum.ent2path.list <- lapply(hum.path.genes.ent,function(gene)hum.allpath.name[str_detect(hum.allpath$entrez_gene_ids,gene)])
names(hum.ent2path.list) <- hum.path.genes.ent
save(hum.allpath.list,file = 'Hum_allpath_list.Rdata')
save(hum.ent2path.list,file = 'Hum_Ent2path_list.Rdata')

mus.allpath <- fread('CPDB_pathways_mouse_genes.tab')
mus.allpath.name <- paste(mus.allpath$pathway,mus.allpath$external_id,mus.allpath$source,sep = '_')
mus.allpath.list <- lapply(mus.allpath$entrez_gene_ids,str2list)
names(mus.allpath.list) <- mus.allpath.name
mus.path.genes.ent <- unlist(mus.allpath.list) %>% unique()
mus.ent2path.list <- lapply(mus.path.genes.ent,function(gene)mus.allpath.name[str_detect(mus.allpath$entrez_gene_ids,gene)])
names(mus.ent2path.list) <- mus.path.genes.ent
save(mus.allpath.list,file = 'Mus_allpath_list.Rdata')
save(mus.ent2path.list,file = 'Mus_Ent2path_list.Rdata')
# all ppi information from pathway commons----------------------------------------------
hum.allppi <- fread('PathwayCommons.8.All.BINARY_SIF.hgnc.txt.sif',header = F)
hum.allppi <- hum.allppi[V1 %in% hum.info$Symbol.hum][V3 %in% hum.info$Symbol.hum]
hum.allppi <- merge(hum.allppi, hum.info[,.(Symbol.hum,Accession.hum)], 
                    by.x = 'V1', by.y = 'Symbol.hum',all.x = T, allow.cartesian = T) %>%
  merge(hum.info[,.(Symbol.hum,Accession.hum)], 
        by.x = 'V3', by.y = 'Symbol.hum',all.x = T, allow.cartesian = T )
hum.allppi <- hum.allppi[,.(Accession.hum.x,Accession.hum.y,V2)] %>% setnames('V2','Interaction_type')
setorder(hum.allppi,na.last = T)
hum.allppi <- unique(hum.allppi,by = c('Accession.hum.x','Accession.hum.y'))
save(hum.allppi,file = 'hum.allppi.Rdata')

hum.allppi.reverse <- data.frame(Accession.hum.x = hum.allppi$Accession.hum.y,
                                 Accession.hum.y = hum.allppi$Accession.hum.x,
                                 # Interaction_type = hum.allppi$Interaction_type,
                                 stringsAsFactors = F)
hum.allppi.dual <- inner_join(hum.allppi,hum.allppi.reverse)[,c('Accession.hum.x','Accession.hum.y')]
hum.allppi.uni <- anti_join(hum.allppi,hum.allppi.reverse)[,c('Accession.hum.x','Accession.hum.y')]

for(i in seq(nrow(hum.allppi.dual))){
  hum.allppi.dual[i,] <- sort(as.character(hum.allppi.dual[i,]))
}
hum.allppi.dual <- filter(hum.allppi.dual,!duplicated(hum.allppi.dual))
hum.allppi.uni <- rbind(hum.allppi.uni,hum.allppi.dual) %>% as.data.table() %>% 
  merge(hum.allppi,by = c('Accession.hum.x','Accession.hum.y'))

save(hum.allppi.uni,file = 'Hum_allppi_uni.Rdata')



# Go information ----------------------------------------------------------
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(topGO)
go.term <- Term(GOTERM)
goterm <- data.table(ID = names(go.term), Name = unname(go.term))
save(goterm,file = 'GOterm.Rdata')

hum.ent2go.list <- as.list(org.Hs.egGO) %>% lapply(names)
mus.ent2go.list <- as.list(org.Mm.egGO) %>% lapply(names)
hum.go2ent.list <- as.list(org.Hs.egGO2EG) %>% lapply(unname)
mus.go2ent.list <- as.list(org.Mm.egGO2EG) %>% lapply(unname)

save(hum.ent2go.list, file = 'Hum_Ent2GO_list.Rdata')
save(hum.go2ent.list, file = 'Hum_GO2Ent_list.Rdata')
save(mus.ent2go.list, file = 'Mus_Ent2GO_list.Rdata')
save(mus.go2ent.list, file = 'Mus_GO2Ent_list.Rdata')

