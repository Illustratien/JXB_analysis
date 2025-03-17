# create correlation network to validate the hypothesis
rm(list=ls())
pacman::p_load(dplyr,purrr,ggplot2,magrittr,foreach)
options(dplyr.summarise.inform = FALSE)

source('src/fun/sssp_anova.R')
# batch test -------------------------------------------------------------------------
merge_dat <- read.csv("data/merge_profile_rep.csv") %>% 
  select(timeid:appl,g:trt,Trait,rep,namCombine) %>% 
  filter(namCombine%in%c("TT61","TT87","tiller_61","spike_number_87")) %>% 
  tidyr::pivot_wider(names_from=namCombine,values_from=Trait) %>% 
  mutate(T61_T87=TT61/TT87,
         SN_T61=spike_number_87/tiller_61) %>% 
  select(-c(tiller_61:TT87)) %>%
  tidyr::pivot_longer(names_to="namCombine",values_to="Trait",cols=c(T61_T87,SN_T61)) %>% 
  bind_rows(.,read.csv("data/merge_profile_rep.csv")) %>% 
  mutate(rep=factor(rep))

data_list <- merge_dat %>%
  filter(!grepl("(L(A)?I$|m(ean|ax)|GBD|GCD87|RUE|dd|5(5|9)|biomasscut)",namCombine)) %>% # the number is not distributed evenly for each date
  group_by(DFG_year,namCombine) %>%
  group_split() %>%
  .[map_lgl(.,~{(nrow(.x)==96&.x$DFG_year[1]=='DFG2019')|
      (nrow(.x)==192&(!.x$DFG_year[1]=='DFG2019'))})] %>%
  # exclude the avereaged trait
  .[map_lgl(.,~{!any(is.na(.x$rep))})] 

pb = txtProgressBar(min = 0, max = length(data_list),
                    style = 3,    # Progress bar style (also available style = 1 and style = 2)
                    width = 30,initial = 0)

sssp_table_ls <- data_list %>%
  imap_dfr(.,~{
    setTxtProgressBar(pb,.y)
    dff <- tryCatch(
      {aov_tb(.x)},
      error=function(cond) {
        print(.y)
        print(cond)
      }
    )
    setTxtProgressBar(pb,.y)
    return(dff)
  })

# -------------------------------------------------------------------------
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

system.time(
  res <- foreach(
    i  = 1:length(data_list),
    df = data_list,
    .packages = c('dplyr','agricolae','tibble','purrr')
  ) %dopar% {
    dff <- tryCatch(
      {hsd_tb(df)},
      error=function(cond) {
        print(i)
        print(cond)
        return(NULL)
      })
    return(dff)
  }
)

invisible(gc())
parallel::stopCluster(my.cluster)
invisible(gc())


resdf <- res %>% 
  .[!map_lgl(.,~{is.null(.x)})] %>% 
  map_dfr(.,~{.x})

write.csv(resdf,"result/HSDtable.csv",row.names = F)
