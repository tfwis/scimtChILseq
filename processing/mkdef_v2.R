
options(warn = 0) 
suppressPackageStartupMessages(library(tidyverse))

indexfile <- commandArgs(trailingOnly=TRUE)[1]
wdir <- commandArgs(trailingOnly=TRUE)[2]
run <- commandArgs(trailingOnly=TRUE)[3]

if(is.na(wdir)) wdir <- getwd()
if(is.na(run)) run <- sub("_.*","",basename(wdir))

df <- readxl::read_xlsx(indexfile,sheet = 1,col_names = TRUE,range = "A3:H387") %>%
  select(1,2,3,5,7) %>%
  rename(mix = 2,lib_name = 3,i5 = 4,i7 = 5) %>%
  fill(lane,.direction = "down") %>%
  fill(mix,.direction = "down") %>%
  filter(
    !is.na(lib_name)
    ) %>%
  mutate(
    mix = as.integer(factor(mix,levels = unique(mix))),
    lane = ifelse(lane == 1,"L001","L002"),
    lib_id = sub(".*-","",lib_name)
    ) %>%
  unite(i5i7,i5,i7,remove = TRUE,sep = "") %>%
  split(.$mix)

for(i in names(df)){
  i_ <- as.integer(i)
  t5 <- readxl::read_xlsx(indexfile,sheet = i_+1,col_names = TRUE)
  if(ncol(t5)==4){
    t5 <- t5 %>% select(2,4) %>% 
      rename(t5_a = 1,t5_b = 2) %>%
      mutate(multi = TRUE)
  }else{
    t5 <- t5 %>% select(2) %>% 
      rename(t5_a = 1) %>%
      mutate(t5_b = NA,multi = FALSE)
  }
  df[[i]] <- expand_grid(df[[i]],t5)
}
cbs <- bind_rows(df) %>%
  mutate(
    t5 = t5_a,
    run = sprintf("R%s",run)
    ) %>%
  unite(cb_a,t5_a,i5i7,remove = FALSE) %>%
  unite(cb_b,t5_b,i5i7,remove = FALSE) %>%
  mutate(cb_b = ifelse(multi,cb_b,NA)) %>%
  select(-t5_a,-t5_b,-lib_name) %>%
  unite(cell,run,lane,mix,lib_id,t5,remove = F) %>%
  split(.$lane)

for(i in names(cbs)){
  la <- sub("L00","lane",i)
  out <- sprintf("%s/%s/def_cellbc.tsv",wdir,la)
  write_tsv(cbs[[i]],out)
}

