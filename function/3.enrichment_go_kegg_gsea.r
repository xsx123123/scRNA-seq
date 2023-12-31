# Function1: Create dir
create_dir <- function(i){
  if(! dir.exists(i)){
    dir.create(i,recursive=T)}
  print_color_note("create dir finish 🍻🍻🍻🍻🍻!!!!")
}
# Function2: ID convert
id_convert <- function(gene_list){
  # gene id convert
  # SYMBOL id convert ENTREZID
  gene.kegg.entrez_id = mapIds(x = org.Hs.eg.db,keys = gene_list,keytype = "SYMBOL",column = "UNIPROT")
  # remove NA  gene
  gene.kegg.entrez_id = na.omit(gene.kegg.entrez_id)
  # https://www.genome.jp/kegg/catalog/org_list.html
  change_gene.kegg_id <- bitr_kegg(gene.kegg.entrez_id,				 
                                   fromType = "uniprot",
                                   toType = 'kegg',
                                   organism = KEGG_database)
  return(change_gene.kegg_id)
}
# Function3: Kegg function
kegg_gene_list <- function(change_gene.kegg_id){
  # change_gene.kegg_id <- kegg_id$kegg 
  change_gene.KEGG <- enrichKEGG(gene = change_gene.kegg_id,
                                 keyType = "kegg",
                                 organism = "hsa",
                                 pvalueCutoff = KEGG_PVALUE,
                                 pAdjustMethod = "BH")
  return(change_gene.KEGG)
}
# Function4: Define kegg description deal
deal_description <- function(ek.rt){
  list <- c()
  for (i in ek.rt$Description){
    new_i <- strsplit(i, " - ")[[1]][1]
    list <- append(list,new_i)}
  return(list)
}
# Function5: Draw kegg plot
kegg <- function(change_gene.KEGG,name,save_folder){
  ek.rt <- as.data.frame(change_gene.KEGG)
  ek.rt$Description <- deal_description(ek.rt)
  ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
  ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
  ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
  ek.rt20 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 20)
  ek.rt20_sort = arrange(ek.rt20, desc(Count))
  temp_data <- ek.rt20_sort %>% mutate(Description = fct_reorder(Description, Count ))
  ek.rt20_1 <- data.frame()
  for (i in c(1:dim(ek.rt20)[1])){
    temp <- as.data.frame(ek.rt20[i,])
    if (is.na(temp$qvalue)){
      temp$qvalue =1
    }else{
      temp$qvalue =temp$qvalue
    }
    ek.rt20_1 <- rbind(ek.rt20_1,temp)
  }
  # bot plot
  p1 <- ggplot(data = ek.rt20_1,aes(x=enrichment_factor,y=Description)) + 
    geom_point(aes(size=Count,color=qvalue)) +
    scale_color_gradient(low="#fc4a1a",high = "#f7b733") + 
    labs(color="pvalue",size="Count",x="Enrichment Factor",y="KEGG term",title="The Most Enriched KEGG Terms") + 
    theme_test() +
    theme(plot.title = element_text(size = 7, face = "bold",hjust = 0.5),
          axis.title.x = element_text(size = 6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text=element_text(size = 6,face = "bold", color="gray50"),
          legend.key.size = unit(0.35, "cm"),
          legend.title = element_text(size = 3, face = "bold",hjust = 0.5),
          legend.text = element_text(size = 3, color = "gray50", face = "bold"))
  p2 <- ggplot(data=temp_data,aes(x=Description,y=Count,fill=enrichment_factor)) +
    geom_bar(stat="identity",width=0.8) + 
    scale_fill_gradient(low="#ffc0cb",high = "#800080") +
    coord_flip() + 
    theme_test() + 
    xlab("KEGG term") + 
    ylab("Gene count") +
    labs(title = "The Most Enriched KEGG Terms") +
    theme(plot.title = element_text(size = 7, face = "bold",hjust = 0.5),
          axis.title.x = element_text(size = 6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size = 6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text=element_text(size = 6,face = "bold", color="gray50"),
          legend.key.size = unit(0.35, "cm"),
          legend.title = element_text(size = 3, face = "bold",hjust = 0.5),
          legend.text = element_text(size = 3, color = "gray50", face = "bold"))
  p <- ggarrange(p1,p2,ncol=2,nrow=1,labels=c("A","B"))
  ggsave(paste0("kegg_",name,".pdf"),plot = p,width = 25,height = 15,units="cm",device="pdf",path=save_folder)
  ggsave(paste0("kegg_",name,".png"),plot = p,width = 25,height = 15,units="cm",device="png",path=save_folder,dpi=1000)
}
# Function6: KEGG analysis
KEGG <- function(gene_list,save_folder,name){
  # test parameter
  # gene_list <- pvalue_up_gene_list
  # name <- "pvalue_down_kegg"
  # save_folder <- save_folder_KEGG_DOWN
  # KEGG analysis
  kegg_id <- id_convert(gene_list)
  kegg_result <- kegg_gene_list(kegg_id$kegg)
  # if KEGG  result == NULL
  if(is.null(kegg_result)){
    cat("\n")
    cat(bgRed("WARRING!!!!!","\n"))
    cat(bgRed("Enrichment KEGG don't enrichment eveyone pathway","\n"))
    cat("\n")
  }else{
    kegg_result_id <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write.csv(kegg_result_id,file.path(save_folder,paste0(name,"-KEGG-ALL-RESULT.csv")))
    kegg_result_id_1 <- as.data.frame(kegg_result_id)
    save_signle_KEGG_data(kegg_result_id_1,save_folder)
    kegg(kegg_result,name,save_folder)
  }
}
# Function7: Get GSEA infor 
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
# Function8: GetFromNamespace
gseaScores <- getFromNamespace("gseaScores", "DOSE")
# Function9: Gseaplot2_fix
gseaplot2_fix <- function(x, geneSetID, title = "", color="green", base_size = 11,
                          rel_heights=c(1.5, .5, 1), subplots = 1:3,pvalue_table = FALSE, ES_geom="line") {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL ## to satisfy codetool
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                          size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                           size=1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer +
    theme(legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  if (length(geneSetID) == 1) {
    ## geneList <- gsdata$geneList
    ## j <- which.min(abs(geneList))
    ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
    ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]
    ## v <- sort(c(v1, v2))
    ## inv <- findInterval(geneList, v)
    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,xmax=~xmax,ymin=~ymin,ymax=~ymax,fill=~I(col)),data=d,alpha=.9,inherit.aes=FALSE)
  }
  ## p2 <- p2 +
  ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
  ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
  ## theme(legend.position="none") +
  ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),color="grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +xlab("Rank in Ordered Dataset") +
    theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  # fix content star
  # get annotation contenet
  anno <- x[geneSetID, c("NES", "pvalue", "p.adjust")]
  colnames(anno) <- c("NES","P-value","P.adjust")
  lab <- paste0(names(anno), "=",  round(anno, 3), collapse=" ")
  # add contenet for p.posplot
  p.pos <- p.pos + annotate("text", x=max(df2$x)/2, y=max(df2$y)*0.8, label= lab,size=2)
  # fix content end
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    rownames(pd) <- pd$Description
    pd <- pd[,-1]
    pd <- round(pd, 4)
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  p.res <- p.res + theme(plot.title = element_text(hjust = 0.5,size=9,face = "bold"))
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())
  if (length(subplots) == 1)
    return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                    l=.2, unit="cm")))
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  p_all <- plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
}
# Function10: Draw anyone pathway GSEA  plot
draw_signle_GSEA_plot <- function(GSEA_KEGG){
  GSEA_KEGG_DATA <- as.data.frame(GSEA_KEGG)
  # GSEA KEGG data deal
  GSEA_KEGG_DATA$Description <- deal_description(GSEA_KEGG_DATA)
  for(i in c(1:dim(GSEA_KEGG_DATA)[1])){
    name <- GSEA_KEGG_DATA$Description[i]
    p1 <- gseaplot2_fix(GSEA_KEGG,i,title= name,base_size = 7)
    ggsave(paste0(name,"_GSEA_enrichment.pdf"),device = "pdf",plot = p1,width=10,height=7,units = "cm",path = save_gsea_folder)
    ggsave(paste0(name,"_GSEA_enrichment.png"),device = "png",plot = p1,width=10,height=7,units = "cm",path = save_gsea_folder,dpi = 1000)
  }
}
# Function11: Save_signle_GSEA_data
save_signle_GSEA_data <- function(GSEA_KEGG_convert_id){
  GSEA_KEGG_convert_id_temp <- as.data.frame(GSEA_KEGG_convert_id)
  GSEA_KEGG_convert_id_temp$Description <- deal_description(GSEA_KEGG_convert_id_temp)
  for(i in c(1:dim(GSEA_KEGG_convert_id_temp)[1])){
    name <- GSEA_KEGG_convert_id_temp$Description[i]
    temp <- GSEA_KEGG_convert_id_temp[i,]
    temp_1 <- t(temp[,-dim(temp)[2]])
    list <- str_split_fixed(temp$core_enrichment, '/', length(strsplit(temp$core_enrichment,"/")[[1]]))
    temp_list <- t(as.data.frame(list))
    rownames(temp_list) <- c(1:length(strsplit(temp$core_enrichment,"/")[[1]]))
    finish_data <- rbind(temp_1,temp_list)
    write.csv(finish_data,file.path(save_gsea_folder,paste0(name,"_KEGG_GASE.csv")))
  }}
# Function12: Save_signle_KEGG_data
save_signle_KEGG_data <- function(kegg_result_id,save_folder){
  kegg_result_id$Description <- deal_description(kegg_result_id)
  for(i in c(1:dim(kegg_result_id)[1])){
    name <- kegg_result_id$Description[i]
    if (grepl("/", name, fixed=TRUE)) {
      name <- gsub("/", "-", name) 
    }else{
      name <- name
    }
    temp <- kegg_result_id[i,]
    temp_1 <- t(temp[,-dim(temp)[2]])
    list <- str_split_fixed(temp$geneID, '/', length(strsplit(temp$geneID,"/")[[1]]))
    temp_list <- t(as.data.frame(list))
    rownames(temp_list) <- c(1:length(strsplit(temp$geneID,"/")[[1]]))
    finish_data <- rbind(temp_1,temp_list)
    write.csv(finish_data,file.path(save_folder,paste0(name,"_KEGG.csv")))
  }}
# Function13: Draw  GO analysis plot
draw_GO_figure <- function(data,name_analysis_type,save_name,save_dir_GO){
  ek.rt <- as.data.frame(data)
  ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
  ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
  ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
  ek.rt20 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 20)
  ek.rt20_sort = arrange(ek.rt20, desc(Count))
  temp_data <- ek.rt20_sort %>% mutate(Description = fct_reorder(Description, Count))
  # bot plot
  p1 <- ggplot(data = ek.rt20,aes(x=enrichment_factor,y=Description)) + 
    geom_point(aes(size=Count,color=qvalue)) +
    scale_color_gradient(low="#fc4a1a",high = "#791E94",
                         breaks=signif(seq(signif(min(ek.rt20$qvalue),3),signif((max(ek.rt20$qvalue)*1.2),3),length.out=5),3)) + 
    scale_x_continuous(limits=c(0,max(ek.rt20$enrichment_factor)*1.2),n.breaks = 10) +
    labs(color="P-value",size="Count",x="Enrichment Factor",y=paste0("GO Terms of ",name_analysis_type),
         title=paste0("The Most Enriched GO Terms of ",name_analysis_type)) + 
    theme_test() +
    theme(plot.title = element_text(size = 6, face = "bold",hjust = 0.5),
          axis.title.x = element_text(size = 6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size =6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text=element_text(size = 4,face = "bold", color="gray50"),
          legend.key.size = unit(0.35, "cm"),
          legend.title = element_text(face = "bold",hjust = 0.5),
          legend.text = element_text(color = "gray50",size=4,face="italic")) +
    guides(fill=guide_legend(nrow=2,title=NULL,override.aes = list(size=6,alpha=0.5)))
  ggsave(file.path(save_dir_GO,paste0("GO_Terms_of_",save_name,"-",name_analysis_type,".pdf")),device = "pdf",plot =p1,width = 10,height = 11,units = "cm")
  ggsave(file.path(save_dir_GO,paste0("GO_Terms_of_",save_name,"-",name_analysis_type,".png")),device = "png",plot = p1,width = 10,height = 11,units = "cm",dpi = 1000)
}
# Function14: Draw all GO analysis plot
draw_all_GO_figure <- function(data,name_analysis_type,save_name,save_dir_GO){
  ek.rt <- as.data.frame(data)
  ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
  ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
  ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
  ek.rt30 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 30)
  temp_data <- ek.rt30 %>% mutate(Description = fct_reorder(Description, type))
  # temp_data sorted
  temp_data$type_name <- factor(temp_data$type_name,level = c("CC","MF","BP"))
  # bot plot
  list <- c(rep("#629460",10),rep("#e42c64",10),rep("#003366",10))
  
  p1 <- ggplot(data = temp_data,aes(x=type_name,y=Description)) + 
    geom_vline(xintercept="CC",lty=5,col="#003366",lwd=0.3,alpha =0.5) +
    geom_vline(xintercept="MF",lty=5,col="#e42c64",lwd=0.3,alpha =0.5) +
    geom_vline(xintercept="BP",lty=5,col="#629460",lwd=0.3,alpha =0.5) +
    geom_point(aes(size=Count,color=qvalue)) +
    scale_color_gradient(low="#fc4a1a",high = "#791E94",
                         breaks=signif(seq(signif(min(temp_data$qvalue),3),
                                           signif((max(temp_data$qvalue)*1.2),3),
                                           length.out=6),3)) + 
    labs(color="P-value",size="Count",x="Enrichment Factor",y=paste0("GO Terms"),
         title=paste0("The Most Enriched GO Terms ")) + 
    theme_test() +
    theme(plot.title = element_text(size = 6, face = "bold",hjust = 0.5),
          axis.title.x = element_text(size = 6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.title.y = element_text(size =6,face = "bold", vjust = 0.5, hjust = 0.5),
          axis.text.x = element_text(size = 7,face = "bold", vjust = 0.5, hjust = 0.5,color=c("#003366","#e42c64","#629460")),
          axis.text.y = element_text(size = 5,face = "bold.italic",color=list),
          axis.text = element_text(size = 4,face = "bold", color="#282523"),
          legend.key.size = unit(0.35, "cm"),
          legend.title = element_text(face = "bold",hjust = 0.5),
          legend.text = element_text(color = "#282523",size=5,vjust = 0,face="bold")) +
    guides(size=guide_legend(nrow=3,override.aes = list(alpha=1),label.vjust=0.5,label.hjust = -0.5))
  # save all GO figure
  ggsave(file.path(save_dir_GO,paste0("GO_Terms_of_",save_name,"-","all",".pdf")),device = "pdf",plot =p1,width = 14,height = 15,units = "cm")
  ggsave(file.path(save_dir_GO,paste0("GO_Terms_of_",save_name,"-","all",".png")),device = "png",plot = p1,width = 14,height = 15,units = "cm",dpi = 1000)
}
# Function15: deal longest GO Description
deal_longest_GO_Description <- function(data){
  for (i in c(1:length(data))){
    # data <- gene.go.BP_data$Description
    # i <- 373
    speace_number <- str_count(data[i], pattern = " ")
    if (speace_number > 5){
      str_split <- strsplit(data[i]," ")[[1]][round(speace_number/2)]
      str_1 <- strsplit(data[i],str_split)[[1]][1]
      str_2 <- strsplit(data[i],str_split)[[1]][2]
      merge_str <- paste0(str_1,"\n",str_2)
      data[i] <- merge_str
    }}
  return(data)
}
# Function16: draw GO other figure
draw_other_plot_of_GO <- function(GO_DATA,GO_NAME,save_dir_GO){
  # draw GO DAG FIGURE
  pdf(file.path(save_dir_GO,paste0(GO_NAME,"_GO_DAG.pdf")),width = 20,height = 20)
  tryCatch({(plotGOgraph(GO_DATA))}, 
           warning = function(w) {cat("draw GOgraph fail!!!!")}, 
           error = function(e) {cat("draw GOgraph fail!!!!")},
           finally = {cat("draw GOgraph !!!!!","\n")})
  dev.off()
  # draw GO DAG FIGURE igraph style
  tryCatch({p1 <- goplot(GO_DATA, showCategory = 22)}, 
           warning = function(w) {p1 <- "NA"}, 
           error = function(e) {p1 <- "NA"})
  if (!exists("p1")){
    cat("DRAW goplot FAIL!!!","\n")
  }else{
    ggsave(file.path(save_dir_GO,paste0(GO_NAME,"_GO_DAG_igraph.pdf")),plot = p1,width = 25,height = 20,device = "pdf")
    ggsave(file.path(save_dir_GO,paste0(GO_NAME,"_GO_DAG_igraph.png")),plot = p1,width = 25,height = 20,device = "png",dpi=600)
  }
  # draw emapplot  plot figure
  GO_DATA2 <- pairwise_termsim(GO_DATA)
  tryCatch({p2 <- emapplot(GO_DATA2)}, 
           warning = function(w) {p2 <- NA}, 
           error = function(e) {p2 <- NA},
           finally = { })
  if(is.na(p2[1])){
    cat("DRAW emapplot FAIL!!!","\n")
  }else{
    ggsave(file.path(save_dir_GO,paste0(GO_NAME,"_GO_emapplot.pdf")),plot = p2,width = 15,height = 10,device = "pdf")
    ggsave(file.path(save_dir_GO,paste0(GO_NAME,"_GO_emapplot.png")),plot = p2,width = 15,height = 10,device = "png",dpi=600)
  }
}
# Function17: GO analysis
GO_analysis <- function(gene_list,save_name,save_dir_GO){
  # test parameter
  # gene_list <- pvalue_up_gene_list
  # save_name <- "UP"
  # save_dir_GO <- save_folder_GO_UP
  # get GO analysis gene list 
  # gene.go <- as.character(gene_list)
  # library annotation package
  # gene id deal
  # gene id convert
  ##############################################
  print_color_note("GO analysis of BP (Biological Process) DO!!!!!!")
  # enrichGO for BP
  gene.go.BP = enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",
                        ont = "BP",pvalueCutoff = GO_PVALUE,qvalueCutoff = GO_QVALUE)
  
  # gene.go.BP == NA
  if (dim(gene.go.BP)[1] == 0){
    cat("\n")
    cat(bgRed("WARRING!!!!!"))
    cat(bgRed("Enrichment BP don't enrichment eveyone pathway"))
    cat("\n")
    gene.go.BP <- as.data.frame(setReadable(gene.go.BP, OrgDb = org.Hs.eg.db, keyType="SYMBOL"))
  }else{
    # ID convert for BP
    gene.go.BP_data <- as.data.frame(gene.go.BP)
    # gene.go.BP <- as.data.frame(setReadable(gene.go.BP, OrgDb = org.Hs.eg.db, keyType="SYMBOL"))
    write.csv(gene.go.BP_data,file = file.path(save_dir_GO,paste0(save_name,"_GO_BP_enrichment_result.csv")))
    # deal longest GO Description
    gene.go.BP_data$Description <- deal_longest_GO_Description(gene.go.BP_data$Description)
    draw_GO_figure(gene.go.BP_data,"BP",save_name,save_dir_GO)
    # create dir
    save_data_dir_BP <- file.path(save_dir_GO,"top30_BP")
    create_dir(save_data_dir_BP)
    # save GO data for BP
    save_everyone_go_data(gene.go.BP_data,"BP",save_data_dir_BP)
    # draw GO dag plot
    draw_other_plot_of_GO(gene.go.BP,"BP",save_dir_GO)
  }
  print_color_note("GO analysis of BP (Biological Process) DONE!!!!!!")
  ##############################################
  print_color_note("GO analysis of CC (Cell Component) DO!!!!!!")
  # enrichGO for CC
  gene.go.CC = enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",
                        ont = "CC",pvalueCutoff = GO_PVALUE,qvalueCutoff = GO_QVALUE)
  # gene.go.CC == NA
  if (dim(gene.go.CC)[1] == 0){
    cat("\n")
    cat(bgRed("WARRING!!!!!","\n"))
    cat(bgRed("Enrichment CC don't enrichment eveyone pathway","\n"))
    cat("\n")
    gene.go.CC <- as.data.frame(setReadable(gene.go.CC, OrgDb = org.Hs.eg.db, keyType="SYMBOL"))
  }else{
    # ID convert for CC
    gene.go.CC_data <- as.data.frame(gene.go.CC)
    # gene.go.CC <- as.data.frame(setReadable(gene.go.CC, OrgDb = org.Hs.eg.db, keyType="SYMBOL"))
    write.csv(gene.go.CC_data,file = file.path(save_dir_GO,paste0(save_name,"_GO_CC_enrichment_result.csv")))
    # deal longest GO Description
    gene.go.CC_data$Description <- deal_longest_GO_Description(gene.go.CC_data$Description)
    draw_GO_figure(gene.go.CC_data,"CC",save_name,save_dir_GO)
    # create dir
    save_data_dir_CC <- file.path(save_dir_GO,"top30_CC")
    create_dir(save_data_dir_CC)
    # save GO data for CC
    save_everyone_go_data(gene.go.CC_data,"CC",save_data_dir_CC)
    # draw GO dag plot
    draw_other_plot_of_GO(gene.go.CC,"CC",save_dir_GO)
  }
  print_color_note("GO analysis of CC (Cell Component) DONE!!!!!!")
  ##############################################
  print_color_note("GO analysis of MF (Molecular Function) DO!!!!!!")
  # enrichGO for MF
  gene.go.MF = enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pAdjustMethod="fdr",
                        ont = "MF",pvalueCutoff = GO_PVALUE,qvalueCutoff = GO_QVALUE)
  # if gene.go.MF == NA
  if (dim(gene.go.MF)[1] == 0){
    cat("\n")
    cat(bgRed("WARRING!!!!!","\n"))
    cat(bgRed("Enrichment MF don't enrichment eveyone pathway","\n"))
    cat("\n")
    gene.go.MF <- as.data.frame(setReadable(gene.go.MF, OrgDb = org.Hs.eg.db, keyType="SYMBOL"))
  }else{
    # ID convert for MF
    gene.go.MF_data <- as.data.frame(gene.go.MF)
    # gene.go.MF <- as.data.frame(setReadable(gene.go.MF, OrgDb = org.Hs.eg.db, keyType="SYMBOL"))
    write.csv(gene.go.MF_data,file = file.path(save_dir_GO,paste0(save_name,"_GO_MF_enrichment_result.csv")))
    # create dir
    save_data_dir_MF <- file.path(save_dir_GO,"top30_MF")
    create_dir(save_data_dir_MF)
    # save GO data for MF
    save_everyone_go_data(gene.go.MF_data,"MF",save_data_dir_MF)
    # draw GO dag plot
    draw_other_plot_of_GO(gene.go.MF,"MF",save_dir_GO)
    # deal longest GO Description
    gene.go.MF_data$Description <- deal_longest_GO_Description(gene.go.MF_data$Description)
    draw_GO_figure(gene.go.MF_data,"MF",save_name,save_dir_GO)
  }
  print_color_note("GO analysis of MF (Molecular Function) DONE!!!!!!")
  ##############################################
  print_color_note("GO analysis of draw all GO result figure DO!!!!!!")
  # sort and save top10 GO:BP
  gene.go.BP[order(gene.go.BP$qvalue),]
  gene.go.BP_top10 <- gene.go.BP[c(1:10),]
  gene.go.BP_top10$type <- 1
  gene.go.BP_top10$type_name <- "BP"
  # sort and save top10 GO:CC
  gene.go.CC[order(gene.go.CC$qvalue),]
  gene.go.CC_top10 <- gene.go.CC[c(1:10),]
  gene.go.CC_top10$type <- 3
  gene.go.CC_top10$type_name <- "CC"
  # sort and save top10 GO:MF
  gene.go.MF[order(gene.go.MF$qvalue),]
  gene.go.MF_top10 <- gene.go.MF[c(1:10),]
  gene.go.MF_top10$type <- 2
  gene.go.MF_top10$type_name <- "MF"
  # get all data
  all_data <- rbind(gene.go.BP_top10,gene.go.CC_top10,gene.go.MF_top10)
  # remove NA
  all_data <- na.omit(all_data)
  # if all_data==NA
  if (dim(all_data)[1]==0){
    cat("\n")
    cat("NOT ENRICHMENT EVERYONE PATHWAYS")
    cat("\n")
  }else{
    # deal longest GO Description
    all_data$Description <- deal_longest_GO_Description(all_data$Description)
    # draw all GO figure
    draw_all_GO_figure(all_data,"all",save_name,save_dir_GO)
  }
  
  print_color_note("GO analysis of draw all GO result figure DONE!!!!!!")
  ##############################################
}
# Function18: remove GO Description NOT suite condition
remove_GO_Description_NOT_suite_condition <- function(Description_list){
  for(i in c(1:length(Description_list))){
    name <- Description_list[i]
    name_deal <- gsub("/", "-", name)
    Description_list[i] <- name_deal
  }
  return(Description_list)
}
# Function19: GO analysis repoert txt
save_everyone_go_data <- function(result,name_dir,save_dir_GO){
  result$Description <- remove_GO_Description_NOT_suite_condition(result$Description)
  GSEA_go_convert_id_temp <- as.data.frame(result)
  GSEA_go_convert_id_temp = GSEA_go_convert_id_temp[order(GSEA_go_convert_id_temp$qvalue),]
  GSEA_go_convert_id_temp_top30 <- GSEA_go_convert_id_temp[c(1:30),]
  for(i in c(1:dim(GSEA_go_convert_id_temp_top30)[1])){
    name <- GSEA_go_convert_id_temp_top30$Description[i]
    temp <- GSEA_go_convert_id_temp_top30[i,]
    temp_1 <- t(temp[,-dim(temp)[2]])
    list <- str_split_fixed(temp$geneID, '/', length(strsplit(temp$geneID,"/")[[1]]))
    temp_list <- t(as.data.frame(list))
    rownames(temp_list) <- c(1:length(strsplit(temp$geneID,"/")[[1]]))
    finish_data <- rbind(temp_1,temp_list)
    write.csv(finish_data,file.path(save_dir_GO,paste0(name_dir,"-",name,"_GO.csv")))
  }}
# Function20: GO analysis repoert txt
GSEA_KEGG <- function(all_deg_data,save_gsea_folder){
  # id convert
  convert_list <- all_deg_data$gene
  gene <- bitr(convert_list,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
  # remove dup content
  gene <- gene[!duplicated(gene[,2]),]
  # chenge rownames
  rownames(gene) <- gene$ENTREZID
  colnames(all_deg_data)[1] <- "SYMBOL"
  # merge data
  GSEA_data <- left_join(gene,all_deg_data,by='SYMBOL');GSEA_data <- GSEA_data[,-1]
  # sort data
  GSEA_data <- GSEA_data %>% mutate(rank = rank(avg_log2FC,  ties.method = "random")) %>% arrange(desc(avg_log2FC))
  gene_list <- GSEA_data$avg_log2FC
  names(gene_list) <- GSEA_data$ENTREZID
  # GSEA analysis
  GSEA_KEGG <- gseKEGG(gene_list, organism = KEGG_database,pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",nPerm = 1000,by= "DOSE")
  # if GSEA_KEGG == NA
  if (dim(GSEA_KEGG)[1]==0){
    cat("\n")
    cat(bgRed("WARRING!!!!!","\n"))
    cat(bgRed("Enrichment GSEA don't enrichment eveyone pathway","\n"))
    cat("\n")
  }else{
    # convert GSEA  result id
    GSEA_KEGG_convert_id <- setReadable(GSEA_KEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    # create GSEA folder
    save_gsea_folder <- save_gsea_folder
    # draw everyone GSEA  plot
    draw_signle_GSEA_plot(GSEA_KEGG)
    # save everyone GSEA PATHWAY DATA
    save_signle_GSEA_data(GSEA_KEGG_convert_id) 
  }
}
# Function21: cell cluster deg result enrichment
RNRICHMENT_FOR_EVERYONE_CELL_CLUSTER_DEG <- function(){
  folder_list <- list.files(file.path(root_dir,deg_file))
  for (i in folder_list){
    # i <- folder_list[1]
    input <- i
    # PATH6 KEGG & GO analysis of UP gene
    print_color_note_type2("up gene kegg analysis (p_value < 0.5 & log2foldchange >1) !!!!")
    # found DEG up gene list (cutoff:p_value < 0.5 & log2foldchange >1)
    deg_folder_name = file.path(root_dir,deg_file,input)
    pvalue_up_file_name = list.files(deg_folder_name,pattern = "*.UP-DEG-result.csv$")
    pvalue_up_deg_result <- read.csv(file.path(deg_folder_name,pvalue_up_file_name))
    pvalue_up_gene_list <- pvalue_up_deg_result$gene
    # GO analysis
    save_folder_GO_UP <- file.path(root_dir,save_name,i,"pvalue_up_GO")
    create_dir(save_folder_GO_UP)
    GO_analysis(pvalue_up_gene_list,"UP",save_folder_GO_UP)
    # KEGG analysis
    save_folder_KEGG_UP <- file.path(root_dir,save_name,i,"pvalue_up_KEGG")
    create_dir(save_folder_KEGG_UP)
    KEGG(pvalue_up_gene_list,save_folder_KEGG_UP,"pvalue_up_kegg")
    # KEGG & GO analysis finish!!!!!!!
    print_color_note_type2("up gene kegg analysis (pvalue < 0.5 & log2foldchange >1) done !!!!")
    ################################# PATH-7 #################################
    # PATH7 KEGG & GO analysis of down gene
    print_color_note_type2("down gene kegg analysis (pvalue < 0.5 & log2foldchange <-1) !!!!")
    # found DEG up gene list (cutoff:p_value < 0.5 & log2foldchange <-1)
    deg_folder_name = file.path(root_dir,deg_file,input)
    pvalue_DOWN_file_name = list.files(deg_folder_name,pattern = "*.DOWN-DEG-result.csv$")
    pvalue_DOWN_deg_result <- read.csv(file.path(deg_folder_name,pvalue_DOWN_file_name))
    pvalue_DOWN_gene_list <- pvalue_DOWN_deg_result$gene
    # GO analysis
    save_folder_GO_DOWN <- file.path(root_dir,save_name,i,"pvalue_down_GO")
    create_dir(save_folder_GO_DOWN)
    GO_analysis(pvalue_DOWN_gene_list,"DOWN",save_folder_GO_DOWN)
    # KEGG analysis
    save_folder_KEGG_DOWN <- file.path(root_dir,save_name,i,"pvalue_down_KEGG")
    create_dir(save_folder_KEGG_DOWN)
    KEGG(pvalue_DOWN_gene_list,save_folder_KEGG_DOWN,"pvalue_down_kegg")
    # KEGG & GO analysis finish!!!!!!!
    print_color_note_type2("down gene kegg analysis (pvalue < 0.5 & log2foldchange <-1) done !!!!")
    ################################# PATH-8 #################################
    # PATH8 GSEA analysis
    # get deg folder dir
    all_deg_dir <- file.path(root_dir,deg_file,i)
    # get all deg file
    all_deg_file <- list.files(all_deg_dir,pattern = "*. DEG-result.csv$")
    # input all deg file
    all_deg_result <- read.csv(file.path(all_deg_dir,all_deg_file))
    # extert name & logFC
    all_deg_data <- all_deg_result[,c("gene","avg_log2FC")]
    all_deg_data <- na.omit(all_deg_data)
    # sort 
    all_deg_data <- all_deg_data[order(-all_deg_data$avg_log2FC),]
    colnames(all_deg_data) <- c("gene","avg_log2FC")
    # save dir
    save_gsea_folder <- file.path(root_dir,save_name,i,"GSEA")
    create_dir(save_gsea_folder)
    # GSEA
    if(isFALSE(GEEA_ANALYSIS)){
      cat("GSEA analysis don't analysis")
    }else{
      GSEA_KEGG(all_deg_data,save_gsea_folder)
    }
  }
}
