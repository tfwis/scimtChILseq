## Remove PCR Duplicates and IVT Variants

suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(furrr))
suppressMessages(library(tidyverse))

plan(multisession(workers = 6))

# wdir <- getwd() 
wdir <- commandArgs(trailingOnly=TRUE)[1]
cellbcs <- commandArgs(trailingOnly=TRUE)[2]
minCount <- 5

rpath <- paste0(wdir, "/sep")
bpath <- paste0(wdir, "/",cellbcs)
opath <- paste0(wdir, "/grs_rmDup/")
opath2 <- paste0(opath,"uniqueReads/")

cellbc <- read_tsv(bpath,show_col_types=FALSE)
cands <- c(cellbc$cb_a,cellbc$cb_b) %>% na.omit()

dir.create(opath)
dir.create(opath2)

grs_uni <- list()
nRead <- list()
bams <- list.files(rpath,full.names = TRUE,pattern = "bam")
bams <- rev(bams)
for(i in 1:length(bams)){
  gr <- readGAlignmentPairs(bams[[i]],use.names=T)
  ra <- ranges(gr)
  gr <- GenomicAlignments::first(gr)
  gr <- GRanges(seqnames(gr),ranges = ranges(gr))
  gr$cell_bc <- stringr::str_sub(names(gr),-23,-1)
  names(gr) <- NULL

  gr_uni <- gr %>%
    resize(width = 1,fix = 'start') %>%
    as_tibble() %>% distinct() %>%
    with(GenomicRanges::GRanges(seqnames,ranges = IRanges(start,end),cell_bc = cell_bc))

  chr <- sub(".bam","",basename(bams[i]))
  out1 <- paste0(opath2,chr,"_unique.rds")
  saveRDS(gr_uni,out1)

  grs_uni[[i]] <- gr_uni
  nRead[[i]] <- ra %>% as.data.frame() %>% as_tibble() %>% 
    mutate(cell_bc = str_sub(names,-23,-1)) %>% distinct(cell_bc,start,end) %>% 
    count(cell_bc)
}

Gr <- do.call(c,grs_uni)
rm(gr_uni); gc()
nFrag <- as_tibble(Gr) %>% 
  dplyr::count(cell_bc) %>%
  filter(n >= minCount)  
nRead <- bind_rows(nRead) %>%
  reframe(n = sum(n),.by = cell_bc) %>%
  filter(n >= minCount)
write_tsv(nFrag,paste0(opath2,'nFrag.tsv'))
write_tsv(nRead,paste0(opath2,'nRead.tsv'))

usecell <- nFrag$cell_bc
Gr <- Gr[Gr$cell_bc %in% usecell]
Gr_cell <- Gr %>% split(Gr$cell_bc)
cbs <- intersect(names(Gr_cell),usecell)
rm(Gr); gc()

cbs_cands <- intersect(cbs,cands)
cbs_cands_sep <- list(
  t5 = str_sub(cbs_cands,1,6),
  i5 = str_sub(cbs_cands,8,15),
  i7 = str_sub(cbs_cands,16,23)
  )
cbs_error <- setdiff(cbs,cands)
cbs_error_sep <- list(
  t5 = str_sub(cbs_error,1,6),
  i5 = str_sub(cbs_error,8,15),
  i7 = str_sub(cbs_error,16,23)
  )
hamdist <- function(type,id){
  stringdist::stringdist(cbs_cands_sep[[type]][id],cbs_error_sep[[type]],method = 'hamming')
}
diff <- function(type,id){
  cbs_error_sep[[type]] != cbs_cands_sep[[type]][id]
}
hd <- map(1:length(cbs_cands),~{
  t5hd <- diff("t5",.x)
  i5hd <- hamdist("i5",.x)>1
  i7hd <- hamdist("i7",.x)>1
  hdsum <- t5hd + i5hd + i7hd
  which(hdsum==0)
  })
uni <- tibble(query = unlist(hd)) %>% dplyr::count(query) %>% filter(n==1)
ec <- tibble(
    corrected = rep(cbs_cands,sapply(hd,length)),
    error = unlist(hd)
    ) %>%
  filter(error %in% uni$query) %>%
  mutate(cell_bc = cbs_error[error]) %>%
  select(cell_bc,corrected)
write_tsv(ec,paste0(opath2,'errorCorrection.tsv'))

ec_ <- tibble(cell_bc = cbs_cands,corrected = cbs_cands) %>% bind_rows(ec)
nFrag_ec <- nFrag %>%
  inner_join(ec_,by = 'cell_bc') %>%
  reframe(n = sum(n),.by = corrected)
nRead_ec <- nRead %>%
  inner_join(ec_,by = 'cell_bc') %>%
  reframe(n = sum(n),.by = corrected)
write_tsv(nFrag_ec,paste0(opath2,'nFrag_ec.tsv'))
write_tsv(nRead_ec,paste0(opath2,'nRead_ec.tsv'))

ecvs <- with(ec,split(cell_bc,corrected))
Gr_cell_ec <- cbs_cands %>%
  map(~{
    if(!is.null(ecvs[[.x]])){
      a <- Gr_cell[[.x]]
      b <- Gr_cell[names(Gr_cell) %in% ecvs[[.x]]]
      names(b) <- NULL
      b <- do.call(c,b)
      res <- c(a,b)
      res$cell_bc <- .x
    }else{
      res <- Gr_cell[[.x]]
      res$cell_bc <- .x
    }
    return(res)
  }) %>%
  structure(names = cbs_cands)
rm(Gr_cell); gc()

idx <- cellbc %>% 
  select(mix,cb_a,cb_b) %>% 
  gather(key,cb,-mix) %>% 
  na.omit() %>%
  mutate(
    key = sub("cb_","",key),
    mix = sprintf("Mix%s",mix)
    ) %>% 
  unite(label,mix,key) %>% 
  with(split(cb,label))
Gr_cell_ec_m <- idx %>%
  map(~{
    Gr_cell_ec[names(Gr_cell_ec) %in% .x] %>%
      structure(names = NULL) %>% do.call(c,.)
    })

for(i in 1:length(Gr_cell_ec_m)){
  lab <- names(Gr_cell_ec_m)[i]
  out <- paste0(opath,'MergedFragments_cands_',lab,'.rds')
  saveRDS(Gr_cell_ec_m[[i]],out)
}

