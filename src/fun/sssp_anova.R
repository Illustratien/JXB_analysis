# functional-------------------------------------------------------------------------
aov_tb<- function(df){
  # basic anova function 
  model<- with(df,aov(Trait ~ rep*nitrogen*appl*timeid*var))
  # summary
  aovtbl<-summary(model) %>% unclass() %>% data.frame() %>% 
    select(Df:Mean.Sq) %>% mutate(F.value=NA,p.value=NA)
  # get all terms 
  terms <- aovtbl %>% rownames() %>% gsub(' ',"",.)
  # nominator 
  nom <- terms%>% 
    .[!grepl("^rep|Resi",.)]
  # position of nominator 
  nomid <-nom %>%  map_dbl(.,~{grep(paste0("^",.x,"$"),terms)})
  # find position of denominator
  denomid <-nom%>% 
    paste0("rep:",.) %>% 
    map_dbl(.,~{grep(paste0("^",.x,"$"),terms)})
  
  for(x in 1:length(nom)){
    # modify f value 
    f <- aovtbl$Mean.Sq[nomid[x]]/aovtbl$Mean.Sq[denomid[x]]
    # calculate p value based on f value with upper tail 
    p <- pf(f, aovtbl$Df[nomid[x]], aovtbl$Df[denomid[x]], lower.tail = F)
    # assign value 
    aovtbl[nomid[x],"F.value"] <- f
    aovtbl[nomid[x],"p.value"] <- p
  }
  # return the targeted rows
  res <- aovtbl[!is.na(aovtbl$F.value),] %>% 
    tibble::rownames_to_column("treatment") %>% 
    mutate(trait=df$namCombine[1],
           Year=df$DFG_year[1],
           part=df$part[1],
           BBCH=df$BBCH[1],
           group=df$group[1],
           treatment=gsub(" ","",treatment))
  return(res)
}

df_trt_modify <- function(df, vec) {
  # create treatment column bsed on significant terms
  vec <- unlist(strsplit(vec,":"))
  vec_sym <- lapply(vec, as.symbol)
  df %>% 
    mutate(trt := interaction(!!!vec_sym, sep = "_"))
}  

hsd_tb<- function(df){
  sig <- aov_tb(df) %>%  dplyr::filter(p.value<.05)
  if(nrow(sig)>0){
    # basic anova function 
    model<- with(df,aov(Trait ~ rep*nitrogen*appl*timeid*var))
    # summary
    aovtbl<-summary(model) %>% 
      unclass() %>% data.frame() %>% 
      select(Df:Mean.Sq) 
    # get all terms 
    terms <- aovtbl %>% rownames() %>% gsub(' ',"",.)
    
    map_dfr(1:nrow(sig),function(x){
      trt_vec <- sig[x,]$treatment
      denomid <- grep(paste0("^rep:",trt_vec,"$"),terms)
      df2 <- df_trt_modify(df,trt_vec)
      hsd_res<- tryCatch(
        {
          with(df2,agricolae::HSD.test(Trait,
                                       trt,
                                       aovtbl$Df[denomid],
                                       aovtbl$Mean.Sq[denomid])) %>%
            .$groups %>% 
            tibble::rownames_to_column("treatment") %>% 
            mutate(combi=trt_vec)
          
        },
        error=function(cond) {
          print(trt_vec)
          print(cond)
          return(NULL)
        })
      
    }) %>%
      .[!is.null(.)] %>% 
      mutate(trait=df$namCombine[1],
             Year=df$DFG_year[1],
             part=df$part[1],
             BBCH=df$BBCH[1],
             group=df$group[1])
    
  }else{return(NULL)}
  
}
