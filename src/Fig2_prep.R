rm(list=ls())
Geno <- c("Alves","Torrild","Esket","Potenzial",
          "Capone","Pionier","Patras","Apertus")
trait.vec<- c("plot_yield","plot_tkw","plot_grain_number","HI",
              "spike_number_87","plot_grain_per_spike","grain_length","grain_width","post61_DM","mobile61_DM","mobile61_DM_rel")
raw <- read.csv("data/merge_profile.csv") %>% 
  filter(!DFG_year=="DFG2019",namCombine%in%trait.vec) %>% 
  mutate(Env=paste(appl,nitrogen,timeid,sep="_"))
dat.sm <- raw %>% group_by(namCombine) %>% group_split()

for (i in c("var","Env")){
  stable_table <- purrr::map_dfr(dat.sm,function(dat1){
    # print(i)  # print(trait0)
    trait0 <-dat1$namCombine[1]
    lambda <- dat1%>%pull(Trait)%>%quantile(., .95, na.rm=TRUE)
    
    stable_table <- table_stability(dat1, "Trait", i,
                                    c(setdiff(c("var","Env"),i),"DFG_year"),
                                    lambda,normalize = F,unit.correct = T) %>% 
      suppressWarnings()
    # wide to long
    long <- 
      stable_table %>% 
      dplyr::select(Genotype,Mean.Trait,
                    Genotypic.superiority.measure,Ecovalence) %>% 
      tidyr::pivot_longer(.,
                          cols = c(Genotypic.superiority.measure,Ecovalence),
                          names_to="SI_name",
                          values_to="SI_value"
      )%>%
      mutate(SI_name= gsub(pattern = '[:.:]',replacement = ' ',x=SI_name)%>%factor,
             Trait=trait0)%>%
      rename(value=Mean.Trait)
    return(long)
  })%>%mutate(Genotype = factor(Genotype,levels = Geno))
  write.csv(stable_table,sprintf("data/stable_table_%s.csv",i),row.names=F)
}


