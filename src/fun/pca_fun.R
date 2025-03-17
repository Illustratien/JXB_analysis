

yield <- read.csv("data/merge_profile.csv") %>% 
  dplyr::filter(trait=='yield',part=='plot')

g_rank <- toolStability::genotypic_superiority_measure(yield,'Trait',genotype = 'var',
                                                       environment = c('appl','DFG_year','timeid','nitrogen'),
                                                       unit.correct = T) %>% 
  dplyr::select(-Mean.Trait)  %>%  arrange(genotypic.superiority.measure) %>% 
  .$Genotype

m_rank <- toolStability::genotypic_superiority_measure(yield,'Trait',genotype = 'g',
                                                       environment = c('DFG_year','var'),unit.correct = T) %>% 
  
  dplyr::select(-Mean.Trait)  %>%  arrange(genotypic.superiority.measure) %>% 
  .$Genotype

col_geno <- c("#003f5c","#58508d","#8a508f","#bc5090", "#de5a79", "#ff6361","#ff8531", "#ffa600")

# shape_geno <- c(0,1,2,5,6,8,12,13)
# names(shape_geno) <- Geno
shape_year <- c(2,1)
names(shape_year) <- c("wet","dry")



choices = 1:2
scale = 1
pc.biplot = TRUE
obs.scale = 1 - scale
var.scale = scale
var.factor = 1.6# arrrow length
ellipse.prob = 0.68
ellipse.linewidth = .6
ellipse.fill = FALSE
circle = FALSE
circle.prob = 0.68
varname.size = 3
varname.adjust = 1.25
varname.color = "darkgreen"
varname.abbrev = FALSE
axis.title = "PC"
point.size = 2
labels.size = 4
ellipse.type = "confidence"
alpha = .7
var.axes = TRUE

pic_vec <- c("TGW","GN","GCD","YCD","spad","spadd","TT7177","TT717","Tiller","tiller_low")
pic_list<- map(pic_vec,~{
  magick::image_read(sprintf("logo/%s.png",.x)) %>% 
    magick::image_fill(., 'none') %>%
    as.raster(.)
})
names(pic_list) <- pic_vec
ellipse = TRUE
ellipse.alpha = 0.6
labels = NULL

ggbi <-function (pcobj,
                 groups = NULL,  
                 col_shape=NULL,
                 lookup=NULL, 
                 shape_v=NULL,
                 shape_split=FALSE,
                 ...) {
  # ellipse.type =   ellipse.alpha = 
  # groups =  ellipse =
  
  if (length(choices) > 2) {
    warning("choices = ", choices, " is not of length 2. Only the first 2 will be used")
    choices <- choices[1:2]
  }
  svd <- get_SVD(pcobj)
  n <- svd$n
  d <- svd$D
  u <- svd$U
  v <- svd$V
  nobs.factor <- ifelse(inherits(pcobj, "prcomp"), sqrt(n - 1), sqrt(n))
  angle <- circle_chol <- ed <- hjust <-
    mu <- sigma <- varname <- xvar <- yvar <- NULL
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  df.v <- var.factor * df.v
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized ", axis.title, 
                         choices, sep = "")
  }else {
    u.axis.labs <- paste(axis.title, choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%%)", 
                                            100 * d[choices]^2/sum(d^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  } else if (!is.null(lookup)){
    
    df.v$varname <- with(lookup,new[match(rownames(v),ori)])
    # print(df.v$varname)
  }else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust <-  with(df.v, (1 - varname.adjust * sign(xvar))/2)
  if(!is.null(shape_v)){
    df.u$shp <- factor(shape_v)
  }
  
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()+
    ggplot2::theme_test()+
    geom_vline(xintercept=0,color="darkgray")+
    geom_hline(yintercept=0,color="darkgray")+
    theme(axis.ticks = element_line(linewidth  = 0.6, color = "black"),
          panel.border = element_rect(colour = ifelse(is.null(frame), 
                                                      NA, "black"), 
                                      fill = NA, size = 0.6),
          axis.line = element_line(colour = "black", linewidth = 0.6))
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  } else {
    if (!is.null(df.u$groups)) {
      
      if(!is.null(col_shape)){
        if(!is.null(shape_v)){
          g <- g + 
            geom_point(aes(color = groups,shape=shp), 
                       alpha = alpha, stroke = 1.1,
                       size = point.size)+
            scale_color_manual(values=col_shape[[1]])+ 
            scale_shape_manual(values =col_shape[[2]] )+
            guides(shape=guide_legend(title="year"))
          
        }else{
          g <- g + 
            geom_point(aes(color = groups), alpha = alpha, 
                       size = point.size)+
            scale_color_manual(values=col_shape[[1]])
        }
        # scale_shape_manual(values=col_shape[[2]])
      }else{
        g <- g + geom_point(aes(color = groups), alpha = alpha, 
                            size = point.size)
      }
      
    }
    else {
      g <- g + geom_point(alpha = alpha, size = point.size)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    geom <- if (isTRUE(ellipse.fill)) 
      "polygon"
    else "path"
    if (isTRUE(ellipse.fill)) {
      g <- g + stat_ellipse(geom = "polygon", aes(group = groups, 
                                                  color = groups, fill = groups), alpha = ellipse.alpha, 
                            linewidth = ellipse.linewidth, type = "norm", 
                            level = ellipse.prob)
    }
    else {
      g <- g + stat_ellipse(geom = "path", aes(group = groups, 
                                               color = groups), linewidth = ellipse.linewidth, 
                            type = "norm", level = ellipse.prob,show.legend = F)
    }
  }
  
  if (var.axes) {
    arrow_style <- arrow(length = unit(1/2, "picas"), 
                         type = "closed", angle = 15)
    g <- g + 
      geom_segment(data = df.v, aes(x = 0, y = 0, 
                                    xend = xvar, yend = yvar), 
                   arrow = arrow_style, color = varname.color, 
                   linewidth = .2)+ 
      ggrepel::geom_text_repel(
        data = subset(df.v,yvar>=-1&xvar>0) ,
        mapping=aes(label = varname, 
                    x = xvar, y = yvar),seed = 199,
        nudge_x = 2.25-subset(df.v,yvar>=-1&xvar>0)$xvar,
        size=labels.size, parse=T,
        force_pull   = 0, # do not pull toward data points
        direction="y",hjust=0,
        max.iter=1e4,max.time=1,
        color = varname.color)+
      ggrepel::geom_text_repel(
        data = subset(df.v,yvar<(-1)|xvar<0),
        mapping=aes(label = varname, 
                    x = xvar, y = yvar),
        size=labels.size, parse=T,
        direction="x",nudge_x=0.1,seed = 199,
        nudge_y = -2-subset(df.v,yvar<(-1))$yvar,
        color = varname.color)+
      scale_x_continuous(limits = c(-3,4.3))+
      scale_y_continuous(limits = c(-3.8,3.3))+
      geom_text(data=
                  data.frame(xvar=c(4,-2.8,-2.8,4),
                             yvar=c(3.2,3.2,-3.8,-3.8),
                             l=c("I","II","III","IV")),
                aes(xvar,yvar,label=l),color="darkgray",fontface="bold"
                
      )+
      annotation_raster(pic_list$GN, -0.3, 0.3, -2.5, -3.5)+
      annotation_raster(pic_list$TGW, 3, 3.5, 1.8, 2.5)+
      annotation_raster(pic_list$GCD, .7,1.7,-2.5, -3.5)+
      annotation_raster(pic_list$YCD, -.7,-1.7,2.5, 3.5)+
      annotation_raster(pic_list$Tiller, -0.6, -1.6, -2, -3)+
      annotation_raster(pic_list$tiller_low, .5, 1, 2.5, 3.3)+
      annotation_raster(pic_list$spad, 3.6,3.8, -0.5, -1.5)+
      annotation_raster(pic_list$spadd, -2, -2.2, 0.6, 1.5)+
      annotation_raster(pic_list$TT7177, 3.2, 4.2, -0.3, .2)+
      annotation_raster(pic_list$TT717, -2, -2.5, -0.3, .12)
    if(shape_split){
      g <- g+facet_grid(~shp)+theme(strip.background = element_blank())
    }
  }
  return(g)
}
pca.check<- function(tr,shape_split=FALSE,typ="geno"){
  df_data <- read.csv("data/merge_profile.csv") %>% 
    filter(namCombine%in%tr,!DFG_year=="DFG2019") %>% 
    group_by(DFG_year,timeid,appl,nitrogen,var,namCombine) %>% 
    summarise(Trait=mean(Trait)) %>% 
    tidyr::pivot_wider(names_from = 'namCombine',values_from = 'Trait') %>% 
    ungroup() %>% 
    mutate(Env=paste(appl,nitrogen,timeid,sep="_")) %>% 
    relocate(.,Env)
  p.mat <- corrplot::cor.mtest(df_data %>% 
                                 .[,-c(1:6)],
                               method = "pearson")$p
  pc <- df_data %>% 
    .[,-c(1:6)] %>% 
    cor(.) %>% 
    ggcorrplot::ggcorrplot(., 
                           type="lower",lab_col = "white",lab_size=4,
                           p.mat = p.mat, sig.level = 0.05,
                           insig = "pch",tl.cex =10,
                           legend.title = parse(text="italic(r)"),
                           color=c( "#CB7000","white", "#007878"),
                           lab = TRUE)+
    theme(legend.title     = element_text(size=8),
          legend.text      = element_text(size=8),
          legend.key.width = unit(dev.size()[1] / 100, "inches"))+
    scale_x_discrete(labels = scales::parse_format())+
    scale_y_discrete(labels = scales::parse_format())
  # geno across year -----------------------------------------------------------------------
  res.pca <-prcomp(df_data %>% dplyr::select(-c(Env:var)),
                   scale=T)

  
  if(identical(setdiff(tr,lookup$ori),character(0))){
    lokup <- lookup
  }else{
    newv<- setdiff(tr,lookup$ori)
    lokup <- rbind(lookup,
                   rep(newv,2) %>% 
                     matrix(.,ncol=2) %>%
                     data.frame() %>% 
                     'colnames<-' (names(lookup))              
    )
  }
  if(typ=="geno"){
    names(col_geno) <- g_rank
    sg <- factor(df_data$var,levels=g_rank)
  }else{
    names(col_geno) <- m_rank
    sg <- factor(df_data$Env,levels=m_rank)
  }

  col_shape <- list(col_geno,shape_year)
  
  
  p<-ggbi(res.pca, 
          groups =sg,
          col_shape=col_shape,
          lookup=lokup,
          shape_v=df_data %>% mutate(y=ifelse(DFG_year=="DFG2020","wet","dry")) %>% .$y,
          shape_split=shape_split
          
  )
  
  lgd <- cowplot::get_plot_component(p+guides(color=guide_legend(ncol=ifelse(typ=="geno",2,1))),
                                     pattern = 'guide-box-right',
                                     return_all=T)
  
  tbla<- res.pca$rotation[,1:2] %>% 
    as.data.frame() %>%
    mutate_all(round_scale) %>%
    `rownames<-`(., with(lokup,new[match(row.names(.),ori)])) 
  
  tbla <- tbla%>% 
    gridExtra::tableGrob(.,
                         theme=gridExtra::ttheme_minimal(core = list(fg_params=list(cex = .7,
                                                                                    fontface="bold")),
                                                         colhead = list(fg_params=list(cex = .7)),
                                                         rowhead = list(fg_params=list(cex = .7,fontface="bold",parse=T)))) %>% 
    # https://stackoverflow.com/questions/43613320/how-to-add-multi-sub-columns-in-gridextratablegrob
    # gtable_add_grob( .,
    #                  grobs = segmentsGrob( name = 'segment',
    #                                        y1 = unit( 0, 'npc' ),
    #                                        gp = gpar( lty = 1, lwd = 1 ) ),
    #                  t = 0, l = 1, r = ncol( .) )%>% 
    gtable_add_grob( .,
                     grobs = segmentsGrob( name = 'segment',
                                           # from y0 to y1 # if y1=y0 = 0, which means horizontal line
                                           y1 = unit( 0, 'npc' ),
                                           # width and type of line 
                                           gp = gpar( lty = 1, lwd = 1 ) ),
                     t = 1,# horizontal position
                     # vertical from left first column to right last column
                     l = 1, r = ncol( .) )%>% 
    gtable_add_grob( .,
                     grobs = segmentsGrob( name = 'segment',
                                           y1 = unit( 0, 'npc' ),
                                           gp = gpar( lty = 1, lwd = 1 ) ),
                     t = -1
                     # nrow(tbla)+1
                     , l = 1, r = ncol( .) )
  
  pp <- cowplot::plot_grid(
    # pc,
    p+theme(legend.position = "none"),
    cowplot::plot_grid(tbla,lgd,ncol=1),
    nrow = 1,rel_widths = c(1,.3)
    # c( 1, 1,.6)
  )
  print(pp)
  return(pp)
}
