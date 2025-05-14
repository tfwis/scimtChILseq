dcolor <- ggsci::scale_color_d3('category20') 
dfill <- ggsci::scale_fill_d3('category20') 
ccolor <- viridis::scale_color_viridis()
ggColorHue <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}

decoyMerge <- function(A,B)
{
  a <- rownames(A)
  b <- rownames(B)
  uni <- union(a,b)
  d1 <- setdiff(a,b)
  d2 <- setdiff(b,a)
  decoy1 <- matrix(0,ncol=ncol(A),nrow=length(d2))
  #decoy1 <- as(decoy1,'sparseMatrix')
  rownames(decoy1) <- d2
  decoy2 <- matrix(0,ncol=ncol(B),nrow=length(d1))
  #decoy2 <- as(decoy2,'sparseMatrix')
  rownames(decoy2) <- d1
  A_ <- rbind(A,decoy1)
  A_ <- A_[uni,]
  B_ <- rbind(B,decoy2)
  B_ <- B_[uni,]
  AB <- cbind(A_,B_)
  colnames(AB) <- c(colnames(A),colnames(B))
  return(AB)
}

ldecoyMerge <- function(A,relabel=F,labels=NULL)
{
  merged <- A[[1]]
  L <- length(A)
  i <- 2
  while (L>=i) {
    merged <- decoyMerge(merged,A[[i]])
    i <- i+1
  } 
  if(relabel){
    if(is_null(labels)){
      colnames(merged) <- paste0('V',1:ncol(merged))
    }else{
      colnames(merged) <- labels
    }
  }
  return(merged)
}

tib2df <- function(tib,RowNames=1) {
  rn <- pull(tib,all_of(RowNames))
  df <- tib %>%
    dplyr::select(-RowNames) %>%
    as.data.frame()
  rownames(df) <- rn
  return(df)
}

relabel_cols <- function(x,prefix = NULL,suffix = NULL,sep = '_'){
  new <- colnames(x)
  if(!is.null(prefix)){
    new <- paste(prefix,new,sep = sep)
  }
  if(!is.null(suffix)){
    new <- paste(new,suffix,sep = sep)
  }
  colnames(x) <- new
  return(x)
}

seu2tab <- function(x,emb = 'umap',usekey = NULL,into=NULL,comp = c(1,2)){
  Tab <- as_tibble(x@meta.data,rownames = 'cell')
  if(!is.null(into)){
    Tab <- Tab %>%
      separate(cell,into=into,remove=F) 
  }
  if(!is.null(emb)){
    if(emb %in% names(x@reductions)){
      Tab_emb <- Embeddings(x[[emb]])[,comp]
      if(is.null(usekey)) usekey <- stringr::str_to_upper(emb)
      colnames(Tab_emb) <- sprintf("%s%s",rep(usekey,length(comp)),comp)
      Tab_emb <- as_tibble(Tab_emb,rownames = 'cell')
      Tab <- inner_join(Tab_emb,Tab,by = 'cell')
      return(Tab)
    }
  }
  return(Tab)
}

style <- function(x,color,reduction='UMAP',axis=c(1,2),Theme=NULL,
                  axis.text=FALSE,axis.title=FALSE,coord_fix =TRUE,
                  palette=TRUE,title=NULL,density=FALSE,
                  legend=TRUE,scale_color_log10=FALSE,size = NULL,
                  guide_pointSize=2.5,...) {
  axis <- paste0(reduction,axis)
  if(is.null(Theme)) Theme <- theme_bw
  g <- ggplot() + Theme()
  if(!density) {
    if(is.null(size)) size <- .1
    if(length(size) == 2){
      g <- g + geom_point(aes_(x=as.name(axis[1]),y=as.name(axis[2]),color=as.name(color),size = as.name(color)),x,...) +
        scale_radius(range = size)
    }else{
      g <- g + geom_point(aes_(x=as.name(axis[1]),y=as.name(axis[2]),color=as.name(color)),x,size = size,...)
    }
  }else{
    if(is.null(size)) size <- .7
    g <- g + geom_density_2d(aes_(x=as.name(axis[1]),y=as.name(axis[2]),color=as.name(color)),x,size = size,...)
  }
  
  if(coord_fix) {
    g <- g + coord_fixed()
  }
  if(!axis.text) {
    g <- g + theme(axis.text = element_blank(),axis.ticks = element_blank())
  }
  if(!axis.title) {
    g <- g + theme(axis.title = element_blank())
  }
  if(is.null(title)) {
    g <- g # + ggtitle(reduction) 
  }else if(!is.na(title)) {
    g <- g + ggtitle(title) 
  }
  if(!legend) g <- g + theme(legend.position = 'none')
  
  if(scale_type(as.matrix(x[,color]))=='discrete' & !is.null(guide_pointSize)) {
    g <- g + guides(color = guide_legend(override.aes = list(size=guide_pointSize)))
  }
  
  if(is.logical(palette)){
    if(palette) {
      if(scale_type(as.matrix(x[,color]))=='discrete') {
        g <- g + ggsci::scale_color_d3('category20')
      }else{
        if(scale_color_log10){
          g <- g + viridis::scale_color_viridis(trans='log10')
        }else{
          g <- g + viridis::scale_color_viridis()
        }
      }
    }
  } else if (!is.logical(palette)) {
    g <- g + palette
  }
  return(g)
}

Red <- scale_color_gradient(low = '#d3d3d3',high = 'red')
dimplot <- function(seu,feature=NULL,reduction='umap',patchwork_nrow=NULL,patchwork_ncol=NULL,red = FALSE,...)  {
  if(length(feature)==0){
    g <- Embeddings(seu[[reduction]]) %>%
      as_tibble() %>% dplyr::rename(V1=1,V2=2) %>%
      mutate(cluster = seu$seurat_clusters) %>%
      style('cluster',reduction = 'V',title = NA,...)
  }else if(any(colnames(seu@meta.data) %in% feature)){
    a <- Embeddings(seu[[reduction]]) %>% as_tibble(rownames='cell') %>% 
      dplyr::rename(V1=2,V2=3) %>% bind_cols(as_tibble(seu@meta.data))
    if(length(feature) == 1){
      g <- style(a,feature,reduction = 'V',title = feature,...) + labs(color="")
    }else if(length(feature) > 1){
      g <- lapply(feature,function(x){
        style(a,x,reduction = 'V',title = x,...) + labs(color="")
      })
      g <- g %>%
        patchwork::wrap_plots(nrow = patchwork_nrow,ncol = patchwork_ncol)
    }
  }else{
    if(red){
      pal <- Red
    }else{
      pal <- ccolor
    }
    if(length(feature)==1){
      g <- Embeddings(seu[[reduction]]) %>%
        as_tibble() %>% dplyr::rename(V1=1,V2=2) %>%
        mutate(exp = FetchData(object = seu, vars = feature)[,1]) %>%
        style('exp',reduction = 'V',title = NA,palette = F,...) + 
        labs(color=feature) + pal
    }else if(length(feature) > 1){
      a <- as_tibble(FetchData(object = seu, vars = feature),rownames = 'cell') %>%
        gather(key = gene,value = Exp,-cell)  %>%
        mutate(gene = factor(gene,levels=feature))
      b <- Embeddings(seu[[reduction]]) %>% 
        as_tibble(rownames = 'cell') %>% dplyr::rename(V1=2,V2=3)
      tmp <- inner_join(b,a,by = 'cell') %>%
        split(.$gene)
      g <- names(tmp) %>%
        lapply(function(x) {
          style(tmp[[x]],color = 'Exp',reduction = 'V',title = x,palette = FALSE,...) + 
            pal + labs(color="")
        }) %>%
        patchwork::wrap_plots(nrow = patchwork_nrow,ncol = patchwork_ncol)
    }
  }
  return(g)
}

vlnplot <- function(seu,feature,group.by=NULL,jitter=FALSE,jitter_size = .01,
                    patchwork_nrow=NULL,patchwork_ncol=NULL,Layer = 'data',
                    Theme = theme_minimal,nonzero=FALSE,scale = 'width')  {
  if(is.null(group.by)) {
    if("seurat_cluters" %in% colnames(seu@meta.data)) {
      group.by <- "seurat_cluters"
    }else{
      group.by <- "Idents"
    }
  }
  a <- FetchData(object = seu,vars = feature,layer = Layer) %>%
    as_tibble(rownames = 'cell') %>%
    gather(key = gene,value = Exp,-cell)  %>%
    mutate(gene = factor(gene,levels=feature))
  b <- seu2tab(seu,emb = NULL) %>%
    left_join(enframe(Idents(seu),'cell','Idents'),by = 'cell')
  b <- b[,c('cell',group.by)]
  tmp <- inner_join(b,a,by = 'cell') %>% split(.$gene)
  gg <- names(tmp) %>%
    lapply(function(x) {
      if(nonzero) tmp[[x]] <- tmp[[x]] %>% filter(Exp>0)
      g <- tmp[[x]] %>%
        ggplot(aes_(as.name(group.by),as.name('Exp'),fill = as.name(group.by))) +
        geom_violin(scale = scale) + theme_minimal() + 
        theme(legend.position = 'none',
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank()) +
        ggtitle(x) + labs(y = 'Expression level')
      if(jitter){
        g <- g + geom_jitter(size = jitter_size)
      }
      return(g)
    }) %>%
    patchwork::wrap_plots(nrow = patchwork_nrow,ncol = patchwork_ncol)
  return(gg)
}

ol2sm <- function(ol,COL = 'query',ROW = 'subject',Ct = 'n')
{
  
  ol <- ol %>%
    dplyr::rename(query = as.name(COL),
           subject = as.name(ROW),
           N = as.name(Ct))
  lev_col <- unique(ol$query)
  lev_row <- unique(ol$subject)
  ol2 <- ol %>%
    mutate(
      query = factor(query,levels=lev_col),
      subject = factor(subject,levels=lev_row)
    )
  sm <- ol2 %>%
    with(
      Matrix::sparseMatrix(
        j = as.double(query),
        i = as.double(subject),
        x = N
      )
    )
  colnames(sm) <- lev_col
  rownames(sm) <- lev_row
  return(sm)
}

extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  Names <- names(x)
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  names(x) <- Names
  trim(x)
}

reg2gr <- function(x) {
  ra <- gsub("-","_",x)
  ra <- sub("_",":",ra)
  ra <- GRanges(
    seqnames = sub(":.*","",ra),
    ranges = IRanges(
      start = as.double(sub(".*:","",sub("_.*","",ra))),
      end = as.double(sub(".*_","",ra)))
  )
  return(ra)
}
gr2reg <- function(x) with(x,sprintf("%s-%s-%s",seqnames,start,end))

defcols <- function(n, l=65) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=l, c=100)[1:n]
}

freqfrac <- function(Table,ROW = 'seurat_clusters',FILL = 'sample',
                     pal = TRUE,countlab = 'nCell',asp = .8){
  tmp <- Table %>%
    select(all_of(c(FILL,ROW))) %>%
    rename(!!'sample' := FILL,
           !!'row' := ROW) %>%
    count(sample,row)
  ct <- reframe(tmp,n = sum(n),.by = row)
  g <- list(
    ggplot() + 
      geom_bar(aes(row,n,fill = sample,label = n),tmp,stat='identity',position = 'stack') +
      ggrepel::geom_text_repel(aes(row,n,label = n),ct) +
      labs(x = ROW,y = countlab,fill = FILL),
    ggplot(tmp,aes(row,n,fill = sample)) + 
      geom_bar(stat='identity',position = 'fill') +
      labs(x = ROW,y = 'Fraction',fill = FILL)
  ) %>%
    patchwork::wrap_plots(guides = 'collect') &
    theme_minimal() &
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          aspect.ratio = asp)
  if(pal) g <- g & ggsci::scale_fill_d3()
  return(g)
}
compclus <- function(x,orgdb = NULL,org = "mouse",keyType = "SYMBOL",ont = "BP",fun = "enrichGO",simplify = FALSE){
  if(is.null(orgdb)){
    if(is.character(org)){
      orgdb <- switch(
        org,
        "mouse" = org.Mm.eg.db::org.Mm.eg.db,
        "human" = org.Hs.eg.db::org.Hs.eg.db )
    }else{
      stop("org or orgdb")
    }
  }
  
  if(is.list(x)){
    res <- clusterProfiler::compareCluster(x,fun = fun,OrgDb = orgdb,keyType = keyType,ont = ont)
  }else if(is.vector(x)){
    res <- clusterProfiler::enrichGO(x,OrgDb = orgdb,keyType = keyType,ont = ont)
  }
  if(simplify) res <- clusterProfiler::simplify(res)
  return(res)
}

linner_join <- function(List, by, labels = NULL) {
  merged <- Reduce(function(x, y) inner_join(x, y, by = by), List)
  if (!is.null(labels)) colnames(merged) <- labels
  return(merged)
}

silh <- function(Dist,lab,lev,mergin=100) {
  if(!is.character(lab)) lab <- as.character(lab)
  mod <- model.matrix(~lab - 1)
  colnames(mod) <- sub("lab","",colnames(mod))
  mod <- t(t(mod)/colSums(mod))
  if(!is.matrix(Dist)) Dist <- as.matrix(Dist)
  d_ <- Dist %*% mod
  Sil <- tibble(id = rownames(d_),
                label = lab) %>%
    mutate(
      intra = map2(id,label,~{ d_[.x,.y] }),
      inter = map2(id,label,~{ 
        a <- d_[.x,]
        idx <- which(rank(a) == 1)
        if(colnames(mod)[idx] == .y){
          idx <- which(rank(a) == 2)
        }
        a[idx]
      }),
      nearest = map2(id,label,~{ 
        a <- d_[.x,]
        idx <- which(rank(a) == 1)
        colnames(mod)[idx]
      }),
      sil = map2(intra,inter,~{(.y-.x)/max(.x,.y)}),
      label = factor(label,levels = lev)
    ) %>%
    unnest(sil) %>%
    arrange(-as.double(label),sil) %>%
    mutate(rank = 1:nrow(.),
           rank = map2(label,rank,~{
             idx <- which(rev(lev) == .x)
             .y + mergin*(idx - 1)
           })
    ) %>%
    unnest(rank)
  avg_Sil <- Sil %>%
    reframe(avg_sil = mean(sil),
            Min = min(rank),
            Max = max(rank),
            .by =  label)
  g <- ggplot() +
    theme_minimal() + dfill +
    theme(axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    geom_area(aes(rank,sil,fill = label,group = label),Sil) +
    geom_segment(aes(x = Min,xend = Max,y = avg_sil,yend = avg_sil),
                 avg_Sil,linetype = 'dashed') +
    labs(x = 'ID',y = "silhouette coef.",fill = "") +
    coord_flip()
  res <- list(plot = g,df = Sil)
  return(res)
}

olcount <- function(A,B,min_ol=NULL,symbol=FALSE)
{
  if(is.null(names(A))){ names(A) <- 1:length(A) }
  if(is.null(min_ol)){
    ol <- GenomicRanges::findOverlaps(A,B) %>% as_tibble()
  }else{
    ol <- GenomicRanges::findOverlaps(A,B,minoverlap = min_ol) %>% as_tibble()
  }
  
  f <- function(idx,A){
    if(!is.null(A$name)) {
      A_ <- A$name[idx]
    }else if(!is.null(names(A))){
      A_ <- names(A)[idx]
    }else{
      A_ <- gr2reg(A)[idx]
    }
    return(A_)
  }
  if(symbol) B$name <- B$symbol
  ol %>%
    mutate(
      queryHits = f(queryHits,A),
      subjectHits = f(subjectHits,B)
    ) %>%
    #group_by(query,subject) %>%
    #summarize(N=n()) %>%　ungroup()
    count(queryHits,subjectHits)
}

f_predcuv <- function(x,usecomp = 1:2){
  abs_ <- unique(x$ab)
  usecomp <- sprintf("DM_%s",usecomp)
  df <- x %>% select(all_of(usecomp)) %>% as.data.frame()
  
  er <- try(princurve::principal_curve(as.matrix(df)), silent = TRUE)
  while(class(er) == "try-error"){
    distant <- which.max(rowSums(df^2))
    df <- df[-distant,]
    x <- x[-distant,]
    er <- try(princurve::principal_curve(as.matrix(df)), silent = TRUE)
  }
  cuv1 <- princurve::principal_curve(as.matrix(df))
  proj <- princurve::project_to_curve(as.matrix(df),cuv1$s[cuv1$ord,])
  res <- list(cuv1 = cuv1,proj = proj)
  return(res)
}

f_lambdaext <- function(x,res,direction = 'RNAPII'){
  cuv1 <- res$cuv1
  proj <- res$proj
  
  if(nrow(x) != length(proj$lambda)){
    idx <- as.integer(rownames(cuv1$s))
    x <- x[idx,]
  }
  
  x <- x %>%
    mutate(lambda = proj$lambda)
  pol_med <- x %>% 
    filter(ab == direction) %>%
    pull(lambda) %>% median()
  # wm <- which.min(
  #   c(abs(min(x$lambda) - pol_med),
  #     abs(max(x$lambda) - pol_med))
  #   )
  # if(length(wm)==0) return(x)
  # if(wm == 1) {
  #   x <- x %>%
  #     mutate(lambda = -1*lambda,
  #            lambda = lambda-min(lambda))
  # }
  wm <- sum(x$lambda > pol_med) - sum(x$lambda < pol_med)
  if(wm>0) {
    x <- x %>%
      mutate(lambda = -1*lambda,
             lambda = lambda-min(lambda))
  }
  return(x)
}

f_preddens <- function(x){
  abs_ <- unique(x$ab)
  nreg <- count(x,ab)
  use <- nreg %>% filter(n>1) %>% pull(ab)
  x <- x %>% filter(ab %in% use)
  rang <- seq(min(x$lambda),max(x$lambda),len = 100)
  
  Bw <- bw.SJ(x$lambda)
  tab_dens <- x %>%
    select(ab,lambda) %>%
    nest(.by = ab) %>%
    mutate(
      dens = map(data,~with(.x,density(lambda,bw = Bw))),
      # dens = map(data,~with(.x,density(lambda,bw = bw.ucv(lambda)))),
      pred = map(dens,~tibble(lambda = rang,
                              density = with(.x,approx(x,y,rang))$y))
      ) %>%
    select(ab,pred) %>% unnest(pred) %>%
    mutate(density = ifelse(is.na(density),0,density)) %>%
    reframe(lambda = lambda,density = density/sum(density),.by = ab) %>%
    spread(ab,density) 
  
  unused <- setdiff(abs_,use)
  if(length(unused)!=0){
    for(i in unused){
      tab_dens[[i]] <- 0
    }
  }
  tab_dens <- tab_dens %>%
    mutate(id = row_number())
  return(tab_dens)
}

f_dens2adj <- function(tab_dens){
  abs_ <- tab_dens %>%
    select(-lambda,-id) %>%
    colnames() %>% unique()
  g_df <- t(combn(length(abs_),2)) %>%
      as_tibble() %>%
      mutate(V1 = abs_[V1],V2 = abs_[V2]) %>%
      mutate(wsdist = map2_dbl(
        V1,V2,~transport::wasserstein1d(tab_dens[[.x]],tab_dens[[.y]],
                                        wa = tab_dens[['lambda']],
                                        wb = tab_dens[['lambda']])
        )) %>%
    filter(!is.nan(wsdist))
  abs_ <- unique(c(g_df$V1,g_df$V2))
  bind_rows(
    g_df, 
    g_df %>% rename(V2 = 1,V1 = 2)
    ) %>%
    mutate(wsdist = ifelse(wsdist == 0,max(wsdist),wsdist)) %>%
    bind_rows(
      tibble(V1 = abs_,wsdist = 0) %>% mutate(V2 = V1)
    ) %>%
    pivot_wider(names_from = V2,values_from = wsdist) %>%
    tib2df() %>% .[abs_,abs_] %>% as.dist()
}
  
f_mds <- function(adj){
  mds_res <- tryCatch(
    expr = {
      MASS::sammon(adj,k = 1,trace = F) %>% 
        .$points %>% as_tibble(rownames = 'ab')  
      }, 
    error = function(x) { return(NA) }
  )
  if(all(is.na(mds_res))) return(NA)

  tmp <- with(mds_res,structure(V1,names = ab))
  if(!is.na(tmp['RNAP2'])){
    if(tmp['RNAP2'] < mean(tmp)){
      mds_res$V1 <- -mds_res$V1
    }
  }
  return(mds_res)
}

pcurve_fitting <- function(data,usecomp = 1:2){
  curve <- f_predcuv(data,usecomp = usecomp)
  sqrt_sum <- curve[['cuv1']]$dist
  data <- f_lambdaext(data,curve)
  densities <- f_preddens(data)
  ab_dist <- f_dens2adj(densities)
  mds <- f_mds(ab_dist)
  list(data = data,curve = curve,sqrt_sum = sqrt_sum,
       densities = densities,ab_dist = ab_dist, mds = mds)
}

gt_plot <- function(dat,x,comp = c(1,2),nrow = NA){
  tab <- dat %>%
    filter(symbol %in% x) %>% 
    mutate(symbol = factor(symbol,levels = x)) %>%
    select(symbol,data) %>% unnest(data)
  g <- tab %>%
    ggplot(aes_(as.name(paste0("DM_",comp[1])),as.name(paste0("DM_",comp[2])),
                color = as.name("ab"),shape = as.name("region"))) + 
      geom_hline(yintercept = 0,linetype = 'dashed') +
      geom_vline(xintercept = 0,linetype = 'dashed') +
      geom_point() + scale_size_area(max_size = 3) + 
      theme_minimal() + facet_wrap(~symbol,scales = 'free',nrow = nrow) +
      theme(panel.grid.minor = element_blank(),
            axis.text = element_text(size = 5),
            aspect.ratio = .6)
  return(g)
}

draw_pcurve <- function(df,gene,nrow = NULL) {
  df <- df %>% 
    filter(symbol %in% gene) %>%
    mutate(
      symbol = factor(symbol,levels = gene),
      data1 = map2(data,curve,~{
        .x <- .x %>% select(DM_1,DM_2) %>%
          rename(DM_1_st = DM_1,DM_2_st = DM_2)
        bind_cols(.x,.y[['cuv1']]$s)
        }),
      data2 = map(curve,~{
        .x <- .x[['cuv1']]
        as_tibble(.x$s[.x$ord,])
        })
      )
  df0 <- df %>% select(symbol,data) %>% unnest(data)
  df1 <- df %>% select(symbol,data1) %>% unnest(data1)
  df2 <- df %>% select(symbol,data2) %>% unnest(data2)
  ggplot() +
    geom_path(aes(DM_1,DM_2),df2,size = .7,color = 'grey',linetype = 'dashed') +
    geom_point(aes(DM_1,DM_2),df1,color = 'red',size = .5) +
    geom_segment(aes(DM_1,DM_2,xend = DM_1_st,yend = DM_2_st),df1) +
    geom_point(aes(DM_1,DM_2,shape = region,color = ab),df0) +
    theme_minimal() + facet_wrap(~symbol,scales = 'free',nrow = nrow) +
    theme(aspect.ratio = .6,panel.grid.minor = element_blank(),
          axis.text = element_text(size = 5))
}

lambda_dens <- function(df,gene,nrow = NULL,ymax = NA,xclip = NULL,rnap = TRUE){
  df <- df %>% 
    filter(symbol %in% gene) %>%
    mutate(symbol = factor(symbol,levels = gene)) %>%
    select(symbol,data) %>% unnest(data)
  if(!is.null(xclip)){
    df <- df %>%
      mutate(lambda = ifelse(lambda>xclip,xclip,lambda))
  }
  if(!rnap) df <- df %>% filter(ab != 'RNAPII')
  df %>%
    ggplot(aes(lambda,fill = ab)) + 
      geom_density(
        # aes(y = after_stat(scaled)),
        color = 0,alpha = .6
        ) + 
      theme_minimal() +
      facet_wrap(~symbol,scales = 'free',nrow = nrow) + 
      theme(aspect.ratio = .4,panel.grid.minor = element_blank(),
            axis.text = element_text(size = 5)) +
      coord_cartesian(ylim = c(0, ymax))
}

lambda_dens_ridges <- function(df,gene,nrow = NULL,ymax = NA,xclip = NULL,rnap = TRUE,lev = NULL){
  df <- df %>% 
    filter(symbol %in% gene) %>%
    mutate(symbol = factor(symbol,levels = gene)) %>%
    select(symbol,data) %>% unnest(data)
  if(!is.null(lev)){
    df <- df %>%
      mutate(ab = factor(ab,levels = lev))
  }
  if(!is.null(xclip)){
    df <- df %>%
      mutate(lambda = ifelse(lambda>xclip,xclip,lambda))
  }
  if(!rnap) df <- df %>% filter(ab != 'RNAPII')
  df %>%
    ggplot(aes(lambda,ab,fill = ab)) + 
      ggridges::geom_density_ridges(color = 0,alpha = .6) +
      theme_minimal() +
      facet_wrap(~symbol,scales = 'free',nrow = nrow) + 
      theme(aspect.ratio = .4,panel.grid.minor = element_blank(),
            axis.text = element_text(size = 5)) +
      coord_cartesian(ylim = c(0, ymax))
}

aligned_plot <- function(dat,Symbol,denoised,Mds,Filter = F,Region = NULL,Ab = NULL){
  b <- dat %>% filter(symbol == Symbol) %>% 
    select(symbol,data) %>% unnest(data)
  if(!is.null(Ab)) b <- b %>% filter(ab %in% Ab)
  aa <- b %>%
    select(ab,reg) %>% nest(.by = ab) %>%
    mutate(data = map(data,~pull(.x,reg)),
           data = map2(ab,data,~as_tibble(denoised[[.x]][.y,],rownames = 'reg'))) %>%
    unnest(data) %>% unite(reg,ab,reg) %>%
    tib2df() %>% t() %>%# scale() %>% 
    cbind(Mds,.) %>% as_tibble()
  b <- b %>% select(ab,reg,bin)
  tmp <- aa %>%
    #select(-dominant_day) %>%
    gather(key = feat,value = expr,-(1:2)) %>% na.omit()  %>%
    mutate(feat = gsub('\\.','-',feat)) %>%
    separate(feat,into = c("ab","reg"),sep = '_') %>%
    left_join(b,by = c('reg','ab')) 
  # if(!is.null(Filter)|Filter) {
  #   if(!is.double(Filter) & Filter) Filter <- .3
  #   max_ps <- tmp %>%
  #     reframe(maxp = max(expr),.by = c(ab,reg)) %>%
  #     filter(maxp>Filter) %>%
  #     select(-maxp)
  #   tmp <- inner_join(max_ps,tmp,by = c("ab","reg"))
  # }
  if(!is.null(Region)){
    if(Region == 'up') {
        tmp <- tmp %>% filter(bin <= 0)
      }else if(Region == 'gb'){
        tmp <- tmp %>% filter(bin >= 0)
      }
  }
  tmp %>%
    #mutate(expr = ifelse(expr > 2,2,expr)) %>%
    #mutate(expr = ifelse(expr < -2,-2,expr)) %>%
    ggplot(aes(V1,V2,color = expr)) + 
      geom_point(size = .5,alpha = .6) + facet_grid(ab~bin) + coord_fixed() +
      xlim(c(-7,7)) + ylim(c(-7,7)) + 
      theme_void() + 
      theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
            legend.position = 'bottom') +
      ccolor
}

cutoff <- function(x){
  idx <- sort(colSums(x),decreasing = TRUE)
  use <- names(idx)[1:floor(0.95*ncol(x))]
  x[,use]
}

li2df <- function(li) {
  if(is.null(names(li))) names(li) <- 1:length(li)
  tibble(name = names(li),
         data = li)
}