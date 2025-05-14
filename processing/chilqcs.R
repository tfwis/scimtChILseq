## Quality Control

base <- commandArgs(trailingOnly=TRUE)[1]
def <- commandArgs(trailingOnly=TRUE)[2]
setwd(paste0(base,'/'))

## Setup 

options(warn = 0) 
suppressPackageStartupMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(mixtools))
if(!file.exists("fragments/")) dir.create("fragments/")
if(!file.exists("fragments/bed/")) dir.create("fragments/bed/")
if(!file.exists("fragments/rds/")) dir.create("fragments/rds/")
if(!file.exists("fragments/")) dir.create("fragments/")
if(!file.exists("qc/")) dir.create('qc/')

## Load count data

df <- def %>%
  read_tsv(show_col_types = FALSE) %>%
  select(mix,cell,multi,cb_a,cb_b)
nreadec <- sprintf("%s/grs_rmDup/uniqueReads/nRead_ec.tsv",dir) %>%
  read_tsv(show_col_types = FALSE)
nfragec <- sprintf("%s/grs_rmDup/uniqueReads/nFrag_ec.tsv",dir) %>%
  read_tsv(show_col_types = FALSE)
nread2d <- df %>%
  left_join(nreadec,by = c("cb_a"="corrected")) %>% rename(nRead_a = n) %>%
  left_join(nfragec,by = c("cb_a"="corrected")) %>% rename(nFrag_a = n) %>%
  left_join(nreadec,by = c("cb_b"="corrected")) %>% rename(nRead_b = n) %>%
  left_join(nfragec,by = c("cb_b"="corrected")) %>% rename(nFrag_b = n) %>%
  mutate(test = map2_lgl(nRead_a,nRead_b,~all(is.na(.x),is.na(.y)))) %>%
  filter(!test) %>% select(-test)
write_tsv(nread2d,'meta_allcands.tsv')

### Set lower threshold

dnormmix <- function(x)
  fit$lambda[1]*dnorm(x,fit$mu[1],fit$sigma[1]) + fit$lambda[2]*dnorm(x,fit$mu[2],fit$sigma[2])
nfrags <- nread2d %>%
  select(mix,nFrag_a,nFrag_b) %>%
  gather(type,nFrag,-c(1,2)) %>%
  filter(!is.na(nFrag))

set.seed(1)
models <- nfrags %>%
  nest(.by = c(mix,type)) %>%
  mutate(model = map(data,~normalmixEM(log10(.x$nFrag+1)))) %>%
  select(-data)
thrs <- models %>%
  mutate(thr = map(model,~{
    idx <- which.max(.x$mu)
    thr <- 10^(qnorm(0.01,.x$mu[idx],.x$sigma[idx]))-1
    return(thr)
  })) %>%
  select(mix,type,type,thr) %>%
  mutate(type = sub("nFrag","thr",type)) %>%
  unnest(thr) %>% spread(type,thr)
if(! "thr_b" %in% colnames(thrs)) thrs[['thr_b']] <- NA
write_tsv(thrs,'qc/thresholds.tsv',progress = F)

nread2d2 <- nread2d %>%
  left_join(thrs,by = 'mix') %>%
  mutate(
    thr_a = nFrag_a>thr_a,
    thr_b = nFrag_b>thr_b,
    use = map2_lgl(thr_a,thr_b,~all(.x,.y)),
    use = ifelse(multi,use,thr_a),
    use = ifelse(is.na(use),FALSE,use)
    ) %>%
  select(-thr_a,-thr_b)
nread2df <- nread2d2 %>% filter(use)
write_tsv(nread2df,'meta_assigned.tsv')

thrs_ <- thrs %>% 
  gather(type,thr,-mix) %>%
  mutate(type = sub('thr_','',type))
nread2d_ <- nread2d %>%
  select(mix,nFrag_a,nFrag_b) %>%
  gather(type,nFrag,-mix) %>%
  mutate(type = sub('nFrag_','',type))
g_thrs <- ggplot() +
  geom_histogram(aes(nFrag),nread2d_,color = 'grey') +
  geom_vline(aes(xintercept = thr),thrs_,linetype = 'dashed') +
  theme_minimal() +
  theme(aspect.ratio = .65) +
  scale_x_log10() +
  facet_grid(type~mix)
ggsave('qc/thrs.pdf',g_thrs)


## Remove doublets (upper threshold)

###### Gaussian read number from k cells

dknorm <- function(x,k=1,mu=1000,sigma=1000) dnorm(x,k*mu,sqrt(k)*sigma)

###### Poisson read number from k cells

dkpois <- function(x,k=1,mu=1000) dpois(x,k*mu)

###### Posterior probability (discrete k)

dpost <- function(k,x,lambda,mu,sigma) {
  p <- function(k) dknorm(x,k,mu,sigma)*dpois(k,lambda)
  p(k)/sum(p(1:10)) # normalize
}

# bap <- exp(-0.02*(1:96))
# bap <- bap/sum(bap)

doup <- nread2df %>%
  split(.$mix) %>%
  map(~{
    tb <- .x %>%
      mutate(
        t5 = map2_chr(cb_a,cb_b,~paste(sub("_.*","",.x),sub("_.*","",.y),sep = ':')),
        r = ifelse(multi,nFrag_a + nFrag_b,nFrag_a)
        ) 
    nlib <- tb$cell %>% str_sub(1,-8) %>%
      sub(".*_","",.) %>% unique() %>% length()
    re <- tb %>% reframe(nobs = n(),totr = sum(r),.by = t5) %>%
      mutate(lambda = -log((nlib-nobs)/nlib))
    mread <- mean(tb$r); sread <- sd(tb$r)
    dtb <- tb %>%
      inner_join(re,by="t5") %>%
      mutate(purrr::map2_dfr(r,lambda,~{
        p <- dpost(1:10,.x,.y,mread,sread)
        tibble(dprob = 1-p[1],kMAP = which.max(p))
      })) %>%
      arrange(dprob) %>% 
      mutate(usage=1:n()/n()) %>%
      select(t5,nobs,cell,r,dprob,kMAP)
    return(dtb) 
  })
bind_rows(doup,.id = 'mix') %>%
  write_tsv('qc/doubletProbs.tsv',progress = F)

mread_ <- nread2df %>%
  mutate(r = ifelse(multi,nFrag_a + nFrag_b,nFrag_a)) %>%
  pull(r) %>% mean()
g_doup <- doup %>%
  bind_rows(.id = 'mix') %>%
  ggplot(aes(as.factor(kMAP),r/1000,color=nobs)) + 
  theme_minimal() + 
  geom_jitter(size=0.5) + 
  geom_hline(yintercept = 1:5*mread_/1000,lty=3,color="red") +
  scale_color_viridis_c() + 
  labs(x = "MAP estimate of cell number",y = "Total read [K]") +
  facet_wrap(~mix)
ggsave('qc/doubletCBs.pdf',g_doup)

sus <- doup %>% bind_rows() %>% filter(dprob > 0.5) %>% pull(cell)
nread2dff <- nread2df %>% filter(!cell %in% sus)
write_tsv(nread2dff,'meta_singlet.tsv')

ncell <- nread2dff %>%
  reframe(nCell = n(),
          expMax = 25*length(unique(sub(".*_","",cb_a))),
          .by = mix) %>%
  mutate(recover = nCell/expMax)
write_tsv(ncell,'qc/nCells.tsv',progress = F)

g_nfrag2d <- nread2d2 %>%
  mutate(use = ifelse(cell %in% sus,FALSE,use)) %>%
  ggplot(aes(nFrag_a,nFrag_b,color = use)) + 
  theme_minimal() +
  geom_point(size = .1) + #geom_density2d() +
  scale_x_log10() + scale_y_log10() + coord_fixed() + 
  scale_color_manual(values = c(`FALSE` = "#d3d3d3",`TRUE` = "red")) +
  theme(legend.position = "none") +
  facet_wrap(~mix,ncol = 2) +
  labs(x = "nFrag a",y = "nFrag b")
ggsave('qc/usedCBs.pdf',g_nfrag2d)


## Make bed files

grs_path <- list.files("grs_rmDup/",pattern = 'rds',full.names = T) 
grs_name <- grs_path %>%
  basename() %>% sub(".rds","",.) %>%
  sub("MergedFragments_cands_Mix","",.)
grs_keys <- str_split(grs_name,'_')
names(grs_path) <- grs_name
names(grs_keys) <- grs_name
usecell <- nread2dff %>% 
  select(mix,cell,cb_a,cb_b) %>%
  gather(key = type,value = cb,-(1:2)) %>%
  na.omit() %>%
  mutate(type = sub("cb_","",type))
for(gr_name in grs_name){
  Mix <- grs_keys[[gr_name]][1]
  Type <- grs_keys[[gr_name]][2]
  usecell_ <- usecell %>% 
    filter(mix == Mix,type == Type) %>%
    with(structure(cell,names = cb))
  gr <- readRDS(grs_path[gr_name]) %>%
    subset(cell_bc %in% names(usecell_))
  new <- usecell_[gr$cell_bc]
  names(new) <- NULL
  gr$name <- new
  
  gr_name_ <- sub("_b","_Pair",sub("_a","_RNAPII",gr_name))
  out1 <- sprintf("fragments/rds/fragments_%s.rds",gr_name)
  out2 <- sprintf("fragments/bed/fragments_%s.bed",gr_name)
  saveRDS(gr,out1)
  rtracklayer::export(gr,out2,format = 'bed')
}


