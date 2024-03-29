phosphoR <- function(base,
                     del_data,
                     minCiteCount = 3,
                     removeDels = F,
                     geneList,
                     plotName,
                     w=12,
                     h=10,
                     returnData=F,
                     plotResults=T
){

  
  ## this function is to analyse phosphosite data on ubiquitination / phosphorylation of proteins
  ## while also including info form pfam scan on protein domains
  ## and intrinsic disorder prediction by IUPred2
  ## and downloads alpha fold predicted structure files from https://alphafold.ebi.ac.uk/
  
  ## input is a vector of human gene symbols (eg geneList = c("geneA","geneB","geneC"))
  
  ## requires the Ubiquitin and phosphorylation datasets to have been downloaded from phosphosite https://www.phosphosite.org/
  ## or use these that were downloaded in August 2021
  ## https://github.com/d0minicO/phosphoR/blob/main/Phosphorylation_site_dataset.txt
  ## https://github.com/d0minicO/phosphoR/blob/main/Phosphosite_Ubiquitination_site_dataset.txt
  ## these files should be downloaded and placed in the "base" directory
  
  ##############
  ## PACKAGES ##
  ##############
  
  require(tidyverse)
  require(magrittr)
  require(data.table)
  require(pfamscanr)
  require(stringr)
  require(cowplot)
  require(patchwork)
  
  #devtools::install_github("wmm27/idpr")
  require(idpr)
  
  
  
  ###############
  ## FUNCTIONS ##
  ###############
  
  # source custom functions
  devtools::source_url("https://github.com/d0minicO/phosphoR/blob/main/customFunctions.R?raw=TRUE")
  
  
  ###########
  ## THEME ##
  ###########
  
  # Modified from Damien's CapCompare
  science_theme <- theme(
    panel.grid.major = element_blank(),#element_line(size = 0.5, color = "grey"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size=8,hjust = 0),
    text = element_text(size = 8),  axis.line = element_line(color="black", size = .5),
    axis.line.x = element_line(color="black", size = .5),
    axis.line.y = element_line(color="black", size = .5),
    plot.margin = unit(rep(.5,4), "lines"),
    #panel.border=element_blank(),
    strip.background = element_blank()
  )
  
  # colours for IDPR plot
  customColors =c("darkolivegreen3", "darkorchid1", "grey65")
  
  

  
  
  ###############
  ## DATA LOAD ##
  ###############
  
  ## load the deleted IDs data if needed (or skip to save time)
  if(removeDels){
    cat("Loading the database of deleted IDs that should be located here \n", del_data, " \n")
     del_accs = readRDS(del_data)
    # set keys to help join faster
    setkey(del_accs, "ID")
  } else{
    cat("Skipping removing deleted uniprot accessions! \n")
  }
  
  
  
  ## load the phosphosite data
  
  cat("Loading the phosphosite data that should be located in \n", base, " \n")
  
  ub = read_tsv(paste0(base,"Phosphosite_Ubiquitination_site_dataset.txt")) %>%
    mutate(Type="Ub")
  
  
  phos = read_tsv(paste0(base,"Phosphorylation_site_dataset.txt")) %>%
    mutate(Type="Phos")
  
  
  dat = rbind.data.frame(ub,phos)
  
  ## get a column of the position information by keeping just the numbers from the MOD_RSD column
  dat$position = makeChar(str_extract_all(dat$MOD_RSD,"\\(?[0-9,.]+\\)?"))
  
  
  
  ## add a citation count column (includes all the low and high throughput papers)
  # mutate the NA to 0 in the numeric columns only
  dat %<>%
    mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
    mutate(citeCount = LT_LIT + MS_LIT + MS_CST)
  
  
  ## just keep the human sequences
  dat %<>%
    filter(ORGANISM == "human")
  
  ## now find the full length of the protein and locate where the Ub modification is
  
  ## remove the proteins with a non uniprot ID
  dat %<>%
    filter(!grepl("NP_",ACC_ID))
  
  
  ## filter on desired gene list
  dat %<>%
    filter(GENE %in% geneList)
  
  
  ## get the names and IDs
  ## have some genes that have more than one accessions so just need to keep the first
  names =
    dat %>%
    dplyr::select(GENE,ACC_ID) %>%
    group_by(GENE) %>%
    dplyr::slice(1) %>%
    distinct()
  
  
  ids = unique(names$ACC_ID)
  #ids = sample(ids,1000)
  
  ## remove deleted IDs
  
  # make the ids table into a data table to allow quick join to the deleted IDs
  ids %<>% data.table()
  
  colnames(ids) = "ID"
  
  # set the key
  setkey(ids,"ID")
  
  
  ## remove the deleted IDs if we want to
  if(removeDels){
    deleted = makeChar(ids[del_accs, nomatch = 0])
    # remove any deleted IDs
    ids %<>%
      filter(ID %notin% deleted)
  }
  
  
  
  
  ids %<>% makeChar
  
  # remove ones that are not listed as deleted but throw an error
  ### if you receive an error that the table is empty when getting sequences from uniprot, 
  ## then try adding the ID that failed to this list
  
  otherDels = c(
    "AAA58698",
    "Q9BV68-2",
    "Q96PX9",
    "Q9ULD5-3"
  )
  
  ids %<>%
    as_tibble() %>%
    filter(value %notin% otherDels) %>%
    makeChar
  
  
  
  ## keep just the good IDs
  
  names %<>%
    filter(ACC_ID %in% ids)
  
  
  
  
  
  #########################
  ## MAIN DATA CRUNCHING ##
  #########################
  

  
  #####
  # 1 # 
  ##### download ALPHAFOLD DATA
  
  cat("Now downloading alphafold predicted structure pdb files! \n")
  
  ## make a temp folder for the pdb files
  pdb_dir = paste0(base,"pdb_files/")
  dir.create(pdb_dir,showWarnings = F)
  
  
  out.df = data.frame()
  for(i in 1:nrow(names)){
    
    name = as.character(names[i,1])
    id = as.character(names[i,2])
    
    ### download the alpha fold predicted pdb file from EBI
    acc_url = paste0("https://alphafold.ebi.ac.uk/files/AF-",id,"-F1-model_v1.pdb")
    
    download.file(acc_url, destfile = paste0(pdb_dir,name,"_",basename(acc_url)))
    
    
    #read.pdb(paste0(pdb_dir,basename(acc_url)))
    
    # report progress
    cat(round(i*100/length(ids),2),"% \n")
    
  }
  
  cat("Done! \n")
  out.df %<>% as_tibble()
  
  
  
  
  #####
  # 2 # 
  ##### pull the sequences from uniprot
  
  cat("Now downloading sequences from uniprot! \n")
  
  out.df = data.frame()
  for(i in 1:nrow(names)){
    
    name = as.character(names[i,1])
    id = as.character(names[i,2])
    
    ### get the AA sequence from uniprot
    acc_url = paste0("https://www.uniprot.org/uniprot/",id,".fasta")
    temp_seq=paste0(read.csv(url(acc_url))[,1],collapse="")
    
    temp.df = data.frame(GENE = name,
                         ACC_ID = id,
                         length = nchar(temp_seq),
                         seq = temp_seq)
    
    out.df= rbind.data.frame(out.df,temp.df)
    
    # report progress
    cat(round(i*100/length(ids),2),"% \n")
    
  }
  
  cat("Done! \n")
  out.df %<>% as_tibble()
  
  #closeAllConnections()
  
  
  
  
  #####
  # 3 # 
  ##### use pfamscanr to locate the protein domains


  ### add a fasta column (with a header) to allow passing to pfamscan
  
  out.df %<>%
    mutate(fasta = paste0(">",GENE,"_",ACC_ID,"\n",seq))
  
  cat("Now scanning PfamScan for protein domains! \n")
  
  ## perform the pfamscan
  out.df2 = data.frame()
  for(i in 1:nrow(out.df)){
    
    x = pfamscan(fasta_str = out.df$fasta[i],
                 email="dominic.owens@utoronto.ca",
                 maxchecktime = 1000)
    
    
    ## only want to rbind the pfamscan results that are not empty other wise error
    if(length(x)>0){
      temp =
        x %>%
        filter(sig==1) %>%
        dplyr::select(-align) %>%
        dplyr::select(3,2) %>%
        mutate(start = env$from,
               end = env$to) %>%
        dplyr::select(-env) %>%
        mutate(GENE = out.df$GENE[i],
               ACC_ID = out.df$ACC_ID[i],
               length = out.df$length[i],
               seq = out.df$seq[i])
    } else {
      warning("No pfamscan results for ", out.df$GENE[i]," \n")
      temp = data.frame()
    }
    
    out.df2 = rbind.data.frame(out.df2,temp)
    
  }
  
  cat("Done! \n")
  out.df2 %<>% as_tibble()
  
  
  ## dont need this
  # now join to the ub and phos data
  #left_join(dat,out.df2,by="ACC_ID") %>%
  #  select(1,3,5,citeCount,name,start,end,length,seq)
  
  
  
  
  #####
  # 4 # 
  ##### calculate the intrinsic disorder to plot using IDPR
  
  cat("Now calculating intrinsic disorder! \n")
  
  out.df3 = data.frame()
  for(i in 1:nrow(out.df)){
    
    prot=as.character(out.df[i,1])
    id=as.character(out.df[i,2])
    
    ## calcualte disorder from uniprot ID
    temp_dat =
      iupred(id, iupredType = "long", 
             plotResults = F, proteinName = NA) %>%
      mutate(Prot=!!prot)
    
    ## save the data into main output
    out.df3 = rbind.data.frame(temp_dat,out.df3)
    
  }
  
  cat("Done! \n")
  
  
  
  #######################
  ## PLOTTING / EXPORT ##
  #######################
  
  genes = 
    dat %>%
    arrange(desc(GENE)) %>%
    dplyr::select(GENE) %>%
    makeChar %>%
    unique()

  
  if(!plotResults){
  cat("Plotting not requested!")
    } else if (plotResults){
    cat("Analysis done, moving on to plotting... \n")
    
    plotlist = list()
    for(gene in genes){
      
      test1 = 
        out.df %>%
        filter(GENE==gene)
      
      if(nrow(test1)>1){
        warning("Multiple genes matched fix input")
        next
      }
      
      
      test2 = 
        out.df2 %>%
        filter(GENE==gene) %>%
        mutate(start=as.numeric(start),
               end=as.numeric(end),
               length=end-start,
               mid=start+(length/2))
      
      test3 = 
        out.df3 %>%
        filter(Prot==gene)
      
      
      test4 =
        dat %>%
        mutate(position=as.numeric(position)) %>%
        filter(GENE==gene & citeCount > minCiteCount)
      
      
      ### 1 ###
      
      # Plot pfam domains
      domains <- 
        ggplot(test2) +
        geom_hline(yintercept = 1, col = "blue",size=.5) +
        geom_rect(aes(xmin=start, xmax=end, ymin=0.5, ymax=1.5), fill="grey") + 
        scale_colour_manual(values = "blue") + 
        geom_text(aes(label = name, x = mid ,y=1), size=2.5, col="black", fontface=1) +
        coord_cartesian(ylim = c(0.5, 1.5), xlim=c(0,as.numeric(test1$length))) +
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0),breaks=NULL)+
        labs(x=NULL, y="Domains") + 
        science_theme + 
        theme(
          legend.position="None",
          #axis.text.y=element_blank(), 
          axis.line.y=element_blank(),
          #axis.ticks.y=element_blank(),
          #axis.text.x=element_blank(), 
          axis.line.x=element_blank(),
          axis.ticks.x=element_line(size = .1),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          plot.margin = unit(c(0,.1,0,.1), "lines"),
          panel.background = element_blank()
        )
      
      
      ### 2 ###
      
      # Plot phos and Ub mods
      if(nrow(test4)>1){
        ptms =
          ggplot(test4,aes(x=position,y=0,col=Type))+
          geom_point(size=3)+#aes(size=citeCount)+
          scale_size_continuous(range = c(4,5))+
          coord_cartesian(ylim = c(-0.0125, 0.0125), xlim=c(0,as.numeric(test1$length))) +
          scale_x_continuous(expand=c(0,0))+
          scale_y_continuous(expand=c(0,0),breaks=NULL)+
          labs(x=NULL, y="PTMs") + 
          ggtitle(gene)+
          science_theme + 
          theme(
            #legend.position="None",
            #axis.text.y=element_blank(), 
            axis.line.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_blank(), 
            axis.line.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            plot.margin = unit(c(0,.1,0,.1), "lines"),
            panel.background = element_blank()
          )
      } else {
        ptms = NULL
        
        domains =
          domains +
          ggtitle(gene)
      }
      
      
      ### 3 ###
      
      # Plot IDPR
      disorder <- 
        ggplot(test3, aes(x=Position, y=IUPred2)) + 
        geom_line(aes(y=c(1)),alpha=.7)+
        geom_line(aes(y=c(0)),alpha=.7)+
        geom_line(aes(y=c(0.5)),alpha=.3, linetype="dotdash")+
        geom_line(aes(color=IUPred2))+
        scale_color_gradient2(high = customColors[1], 
                              low = customColors[2], mid = customColors[3], 
                              midpoint = 0.5)+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0),breaks=c(0,.5,1))+
        guides(color=FALSE)+
        xlab(NULL)+
        ylab("Intrinsic Disorder")+
        theme_bw()+
        science_theme + 
        theme(panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x=element_blank(), 
              axis.line.x=element_blank(),
              axis.line.y=element_blank(),
              axis.ticks.x=element_blank())
      
      
      ### 4 ###
      
      # combine plots and save into list or just plot alone
      pg <- plot_grid(ptms, domains, disorder, align = "v", axis = "lrtb", nrow=3, rel_heights = c(.5,.5,2))
      
      if(nrow(test4)<1){
        warning("No modifications to plot")
        pg <- plot_grid(domains, disorder, align = "v", axis = "lrtb", nrow=3, rel_heights = c(.5,2))
      }
      
      
      #############
      ## OUTPUTS ##
      #############
      
      ## 1 ##
      
      ## each individual protein map
      
      plotFolder = paste0(base,plotName,"_plot/")
      dir.create(plotFolder, showWarnings = F)
      
      ggsave(pg,
             filename = paste0(plotFolder,gene,".pdf"),
             width=4,
             height=3)
      
      plotlist[[gene]] = pg
      
    }
    
    
    ## 2 ##
    
    ### combination plot
    
    p = wrap_plots(plotlist) +
      plot_layout(guides="collect")
    
    
    ggsave(p,file=paste0(plotFolder,plotName,"_combo.pdf"),
           width=w,
           height=h)
    
    ggsave(p,file=paste0(plotFolder,plotName,"_combo.png"),
           width=w,
           height=h)
    
  }
  
  ## 3 ##
  
  ### return a list of dfs containing the underlying data (if requested)
  
  if(returnData){
    cat("Returning the data frames as a list \n")
    
    out_list =
      list(dat,
           out.df,
           out.df2,
           out.df3)
    
    
    return(out_list)
  } else {
    cat("Dataframe not requested \n")
    return(p)
  }
  
}