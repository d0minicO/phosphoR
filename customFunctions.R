## useful custom functions from various sources

.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, object.size)
  obj.prettysize <- sapply(obj.size, function(r) prettyNum(r, big.mark = ",") )
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size,obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  out <- out[c("Type", "PrettySize", "Rows", "Columns")]
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (head)
    out <- head(out, n)
  out
}


# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}




`%notin%` = function(x,y) !(x %in% y)

makeChar = function(x) (c(x) %>% unlist %>% unname)


standard_error <- function(x) sd(x) / sqrt(length(x))


geom_uperrorbar <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomUperrorbar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}


GeomUperrorbar <- ggproto("GeomUperrorbar", Geom,
                          default_aes = aes(colour = "black", size = 0.5, linetype = 1, width = 0.5,
                                            alpha = NA),
                          
                          draw_key = draw_key_path,
                          required_aes = c("x", "y", "ymax"),
                          
                          setup_data = function(data, params) {
                            data$width <- data$width %||%
                              params$width %||% (resolution(data$x, FALSE) * 0.9)
                            
                            transform(data,
                                      xmin = x - width / 2, xmax = x + width / 2, width = NULL
                            )
                          },
                          draw_panel = function(data, panel_scales, coord, width = NULL) {
                            GeomPath$draw_panel(data.frame(
                              x = as.vector(rbind(data$xmin, data$xmax, NA, data$x,   data$x)),
                              y = as.vector(rbind(data$ymax, data$ymax, NA, data$ymax, data$y)),
                              colour = rep(data$colour, each = 5),
                              alpha = rep(data$alpha, each = 5),
                              size = rep(data$size, each = 5),
                              linetype = rep(data$linetype, each = 5),
                              group = rep(1:(nrow(data)), each = 5),
                              stringsAsFactors = FALSE,
                              row.names = 1:(nrow(data) * 5)
                            ), panel_scales, coord)
                          }
)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}


elbow_plot <- function(mydata,plot_dir,plot_name){
  ## takes a data matrix
  ## a directory to output the plot
  ## and a label to name the plot in the form "my_clustering"
  
  ## generates an elbow plot for up to 15 clusters in the data
  
  
  ## determine number of clusters
  wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                       centers=i)$withinss)
  ## elbow plot
  
  edat = data.frame(wss,num_clust=1:15)
  
  ggplot(edat, aes(num_clust,wss))+
    geom_line()+
    geom_point()+
    ggtitle("Elbow plot", "for cluster number determination")+
    xlab("Number of Clusters")+
    ylab("Within groups sum of squares")+
    scale_x_continuous(breaks=1:15)+
    theme_bw()
  
  
  ggsave(filename=paste0(plot_dir,plot_name,"_elbow_plot.pdf"),
         width=4,
         height=4)
  
}











get_go <- function (universe, de_genes, genome = "hg18", symbol = "geneSymbol"){#,KEGGID){
  
  # this is a custom function from Joe Harman, MRC WIMM, Oxford
  
  require(goseq)
  #require(KEGG.db)
  #xx <- KEGGID
  i <- as.integer(universe %in% de_genes)
  names(i) <- universe
  pwf <- nullp(i, genome, symbol, plot.fit = F)
  GO_List <- list(MF = goseq(pwf, genome, symbol, test.cats = c("GO:MF")), 
                  BP = goseq(pwf, genome, symbol, test.cats = c("GO:BP")), 
                  CC = goseq(pwf, genome, symbol, test.cats = c("GO:CC")))#, 
  #KEGG = goseq(pwf, genome, symbol, test.cats = c("KEGG")))
  GO_List <- lapply(GO_List, function(GO) {
    if ("term" %in% colnames(GO)) {
      return(GO)
    }
    #else {
    #  KEGG <- GO
    #  term <- sapply(KEGG[, 1], function(ID) {
    #    txt <- xx[names(xx) == ID]
    #    txt <- substr(txt, regexpr("=", txt) + 2, nchar(txt))
    #    return(txt)
    #  })
    #  KEGG <- data.frame(KEGG, term)
    #  return(KEGG)
    #}
  })
  return(GO_List)
}


go_plot <- function(all_subs_go,plotName,plot_dir){
  # function takes a list outputted from get_go
  # and a character string (no spaces) to call the plot 
  # and the plotthing output directory
  # (ie whether this is list of all substrates, or only Nterminal Proline substrates)
  
  # loop through each type of GO term and plot each one
  all_terms = names(all_subs_go)
  
  for(i in 1:length(all_terms)){
    
    # get just one of the GO term types as a df
    go_out = 
      all_subs_go[[i]]
    
    # filter the pvals for significantly overrepresented ones
    go_out %<>%
      filter(over_represented_pvalue<.05)
    
    # calculate gene ratio, clip the p values so logging can happen, take just the top 10 
    go_out %<>%
      mutate(geneRatio = (numDEInCat/numInCat)) %>%
      dplyr::select(term, over_represented_pvalue, numDEInCat, numInCat, geneRatio, category) %>%
      arrange(over_represented_pvalue) %>%
      mutate(over_represented_pvalue=ifelse(over_represented_pvalue==0,2e-16,over_represented_pvalue)) %>%
      mutate(minuslogp=-(log10(over_represented_pvalue))) %>%
      dplyr::slice(1:10)
    
    
    # plot name based on what GO type this is
    pname=paste0("GO_", names(all_subs_go)[i])
    
    
    ## do JH style circe GO term plot
    p <- go_out %>%
      ggplot(aes(x = pname, y = term, size = minuslogp, fill = geneRatio)) +
      geom_point(pch=21) +
      theme_classic() +
      labs(title = paste0(pname,"_", plotName), x = NULL, y = NULL)
    
    ggsave(plot=p,
           filename=paste0(plot_dir,pname, plotName,".pdf"),
           width=7,
           height = 3)
    
    
  }
  
  
}


go_genes <- function(GO_list,p_thresh,ensembl,gene_set){
  ## function takes a GO list from get_go function
  ## p value threshold to use
  ## ensembl biomart object for gene / go matching
  ## the de_genes set of the genes used in get_go
  
  out_df = data.frame()
  for (i in 1:length(GO_list)){
    
    
    ## skip the KEGG terms as their cats are not standard GO and cba
    if(names(GO_list[i])=="KEGG"){
      next
    }
    
    
    # get just the top 10 categories that are significant
    temp = 
      GO_list[[i]] %>%
      filter(over_represented_pvalue<p_thresh) %>%
      dplyr::slice(1:10)
    
    ## now get the list of GO categories
    cats = temp$category
    
    
    test = data.frame()
    for(cat in cats){
      
      
      cat("looking for all genes associated with term ", cat, "\n")
      
      # get the total list of genes associated with these go terms
      temp =
        getBM(attributes = c('external_gene_name', 'go_id', 'name_1006'), 
              filters = 'go', 
              values = cat, 
              mart = mart)
      
      
      test = rbind.data.frame(test,temp)
      
    }
    
    
    ## filter those genes that are in our gene set
    test %<>%
      filter(external_gene_name %in% gene_set) %>%
      filter(go_id %in% cats)
    
    ## add an identifier that allows this go type to be identified
    test %<>%
      mutate(Class_GO = paste0("GO_",names(GO_list[i])))
    
    ## spread into long format to get one gene per row
    test %<>%
      group_by(external_gene_name) %>%
      dplyr::select(-go_id) %>%
      mutate(order = seq_along(external_gene_name)) %>%
      spread(key = order, value = name_1006) %>%
      ungroup()
    
    ## will only work if there are max 3 terms associated with each gene
    ## so just keeping the first 3 go terms
    if(ncol(test)>5){
      test %<>%
        dplyr::select(1:5)
    } else if(ncol(test)==4){
      # if only 2 cols present, create a dummy 3rd column
      test %<>%
        mutate(GO_2="")
    } else if(ncol(test)==3){
      # for safety if only 1 present, create 2 dummy columns
      test %<>%
        mutate(GO_2="",
               GO_3="",)
    }
    
    colnames(test)[3:5] = c("GO_1","GO_2","GO_3")
    
    test %<>%
      mutate(GO_2 = as.character(GO_2),
             GO_3 = as.character(GO_3)) %>%
      mutate(GO_2 = if_else(is.na(GO_2),"",GO_2),
             GO_3 = if_else(is.na(GO_3),"",GO_3))
    
    
    ## make into useful columns only
    
    test %<>%
      mutate(across(contains("GO"), ~ as.character(.x))) %>%
      mutate(Combined_GO_terms = if_else(GO_2=="",GO_1,
                                         if_else(GO_3=="",paste0(GO_1,", ",GO_2),
                                                 if_else(GO_3!="",paste0(GO_1,", ",GO_2,", ",GO_3),"wentWrong")))) %>%
      dplyr::select(-starts_with("GO_")) %>%
      distinct()
    
    
    
    
    ## rbind the output df
    out_df = rbind.data.frame(out_df,test)
    
  }
  return(out_df)
}



## custom theme
theme_dom = function(){
  theme_bw()+
    theme(panel.grid=element_blank())
}
