# author:zhang jian
# date:2023.9.5
# version:1.0v
# description: this is scRNA-seq VENN ENRICHMENT
##----------------------------------###---------------------------------------##
# Function1: Create dir
create_dir <- function(i){
  if(! dir.exists(i)){
    dir.create(i)}
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
  change_gene.KEGG <- enrichKEGG(gene = change_gene.kegg_id,
                                 keyType = "kegg",
                                 organism = KEGG_database,
                                 pvalueCutoff = 0.05,
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
  # change_gene.KEGG <- kegg_result
  ek.rt <- as.data.frame(change_gene.KEGG)
  ek.rt$Description <- deal_description(ek.rt)
  ek.rt <- separate(data=ek.rt, col=GeneRatio, into = c("GR1", "GR2"), sep = "/") 
  ek.rt <- separate(data=ek.rt, col=BgRatio, into = c("BR1", "BR2"), sep = "/") 
  ek.rt <- mutate(ek.rt, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) 
  ek.rt20 <- ek.rt %>% filter(row_number() >= 1,row_number() <= 20)
  ek.rt20_sort = arrange(ek.rt20, desc(Count))
  temp_data <- ek.rt20_sort %>% mutate(Description = fct_reorder(Description, Count ))
  # bot plot
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
  ggsave(paste0("kegg_",name,".pdf"),plot = p2,width = 10,height = 10,units="cm",device="pdf",path=save_folder)
  ggsave(paste0("kegg_",name,".png"),plot = p2,width = 10,height = 10,units="cm",device="png",path=save_folder,dpi=1000)
}
# Function6: KEGG analysis
KEGG <- function(gene_list,name,save_folder){
  # test parameter
  # gene_list <- pvalue_up_gene_list
  # save_folder <- save_folder_KEGG
  # name <- "pvalue_up_kegg"
  # KEGG analysis
  kegg_id <- id_convert(gene_list)
  kegg_result <- kegg_gene_list(kegg_id$kegg)
  # if KEGG  result == NULL
  kegg_data <- as.data.frame(kegg_result)
  if(dim(kegg_data)[1]==0){
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
draw_signle_GSEA_plot <- function(GSEA_KEGG,save_gsea_folder){
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
save_signle_GSEA_data <- function(GSEA_KEGG_convert_id,save_gsea_folder){
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
  ek.rt20[is.na(ek.rt20)] <- 1
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
  temp_data[is.na(temp_data)] <- 1
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
    speace_number <- str_count(data[i], pattern = " ")
    if(speace_number > 4){
      str_split <- strsplit(data[i]," ")[[1]][3]
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
  # save_name <- name 
  # save_dir_GO <- save_folder_GO
  # get GO analysis gene list 
  # gene.go <- as.character(gene_list)
  # library annotation package
  # gene id deal
  # gene id convert
  ##############################################
  print_color_note("GO analysis of BP (Biological Process) DO!!!!!!")
  # enrichGO for BP
  gene.go.BP = enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",
                        ont = "BP",pvalueCutoff = 0.01,qvalueCutoff = 0.05)
  
  # gene.go.BP == NA
  gene.go.BP_data <- as.data.frame(gene.go.BP)
  # gene.go.CC == NA
  if (dim(gene.go.BP_data)[1] == 0){
    cat("\n")
    cat(bgRed("WARRING!!!!!","\n"))
    cat(bgRed("Enrichment BP don't enrichment eveyone pathway","\n"))
    cat("\n")
    save_BP_result <- data.frame(ID=NA,Description=NA,GeneRatio=NA,BgRatio=NA,pvalue=NA,p.adjust=NA,qvalue=NA,geneID=NA,Count=NA,type=NA,type_name=NA)
  }else{
    # ID convert for CC
    save_BP_result <- as.data.frame(setReadable(gene.go.BP, OrgDb = org.Hs.eg.db, keyType="ENSEMBL"))
    write.csv(save_BP_result,file = file.path(save_dir_GO,paste0(save_name,"_GO_BP_enrichment_result.csv")))
    # deal longest GO Description
    save_BP_result$Description <- deal_longest_GO_Description(save_BP_result$Description)
    draw_GO_figure(save_BP_result,"BP",save_name,save_dir_GO)
    # create dir
    save_data_dir_BP <- file.path(save_dir_GO,"top30_BP")
    create_dir(save_data_dir_BP)
    # save GO data for CC
    save_everyone_go_data(save_BP_result,"BP",save_data_dir_BP)
    # draw GO dag plot
    draw_other_plot_of_GO(gene.go.BP,"BP",save_dir_GO)
  }
  print_color_note("GO analysis of BP (Biological Process) DONE!!!!!!")
  ##############################################
  print_color_note("GO analysis of CC (Cell Component) DO!!!!!!")
  # enrichGO for CC
  gene.go.CC = enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",
                        ont = "CC",pvalueCutoff = 0.01,qvalueCutoff = 0.05)
  gene.go.CC_data <- as.data.frame(gene.go.CC)
  # gene.go.CC == NA
  if (dim(gene.go.CC_data)[1] == 0){
    cat("\n")
    cat(bgRed("WARRING!!!!!","\n"))
    cat(bgRed("Enrichment CC don't enrichment eveyone pathway","\n"))
    cat("\n")
    save_CC_result <- data.frame(ID=NA,Description=NA,GeneRatio=NA,BgRatio=NA,pvalue=NA,p.adjust=NA,qvalue=NA,geneID=NA,Count=NA,type=NA,type_name=NA)
  }else{
    # ID convert for CC
    save_CC_result <- as.data.frame(setReadable(gene.go.CC, OrgDb = org.Hs.eg.db, keyType="ENSEMBL"))
    write.csv(save_CC_result,file = file.path(save_dir_GO,paste0(save_name,"_GO_CC_enrichment_result.csv")))
    # deal longest GO Description
    save_CC_result$Description <- deal_longest_GO_Description(save_CC_result$Description)
    draw_GO_figure(save_CC_result,"CC",save_name,save_dir_GO)
    # create dir
    save_data_dir_CC <- file.path(save_dir_GO,"top30_CC")
    create_dir(save_data_dir_CC)
    # save GO data for CC
    save_everyone_go_data(save_CC_result,"CC",save_data_dir_CC)
    # draw GO dag plot
    draw_other_plot_of_GO(gene.go.CC,"CC",save_dir_GO)
  }
  print_color_note("GO analysis of CC (Cell Component) DONE!!!!!!")
  ##############################################
  print_color_note("GO analysis of MF (Molecular Function) DO!!!!!!")
  # enrichGO for MF
  gene.go.MF = enrichGO(gene = gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pAdjustMethod="fdr",
                        ont = "MF",pvalueCutoff = 0.01,qvalueCutoff = 0.05)
  # if gene.go.MF == NA
  gene.go.MF_data <- as.data.frame(gene.go.MF)
  if (dim(gene.go.MF_data)[1] == 0){
    cat("\n")
    cat(bgRed("WARRING!!!!!","\n"))
    cat(bgRed("Enrichment MF don't enrichment eveyone pathway","\n"))
    cat("\n")
    save_MF_result <- data.frame(ID=NA,Description=NA,GeneRatio=NA,BgRatio=NA,pvalue=NA,p.adjust=NA,qvalue=NA,geneID=NA,Count=NA,type=NA,type_name=NA)
  }else{
    # ID convert for MF
    save_MF_result <- as.data.frame(setReadable(gene.go.MF, OrgDb = org.Hs.eg.db, keyType="ENSEMBL"))
    write.csv(save_MF_result,file = file.path(save_dir_GO,paste0(save_name,"_GO_MF_enrichment_result.csv")))
    # deal longest GO Description
    save_MF_result$Description <- deal_longest_GO_Description(save_MF_result$Description)
    draw_GO_figure(save_MF_result,"MF",save_name,save_dir_GO)
    # create dir
    save_data_dir_MF <- file.path(save_dir_GO,"top30_MF")
    create_dir(save_data_dir_MF)
    # save GO data for MF
    save_everyone_go_data(save_MF_result,"MF",save_data_dir_MF)
    # draw GO dag plot
    draw_other_plot_of_GO(gene.go.MF,"MF",save_dir_GO)
  }
  print_color_note("GO analysis of MF (Molecular Function) DONE!!!!!!")
  ##############################################
  print_color_note("GO analysis of draw all GO result figure DO!!!!!!")
  # sort and save top10 GO:BP
  save_BP_result[order(save_BP_result$qvalue),]
  save_BP_result_top10 <- save_BP_result[c(1:10),]
  save_BP_result_top10$type <- 1
  save_BP_result_top10$type_name <- "BP"
  # sort and save top10 GO:CC
  save_CC_result[order(save_CC_result$qvalue),]
  save_CC_result_top10 <- save_CC_result[c(1:10),]
  save_CC_result_top10$type <- 3
  save_CC_result_top10$type_name <- "CC"
  # sort and save top10 GO:MF
  save_MF_result[order(save_MF_result$qvalue),]
  save_MF_result_top10 <- save_MF_result[c(1:10),]
  save_MF_result_top10$type <- 2
  save_MF_result_top10$type_name <- "MF"
  # get all data
  all_data <- rbind(save_BP_result_top10,save_CC_result_top10,save_MF_result_top10)
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
# Function20: deal everyone go kegg of venn data
everyone_venn_data_enrichment <- function(){
  for (i in 1:length(venn_list)){
    # i <- 1
    data_path <- file.path(root_dir,venn_list[i])
    name <- strsplit(venn_list[i],".csv")[[1]][1]
    # print run condition
    print_color_note_type3(paste0(name," KEGG-GO DO!!!!"))
    # input veen data
    gene_list <- read.csv(data_path)
    gene_list <- gene_list$gene_symbol
    # GO analysis
    save_folder_GO <- file.path(save_dir,paste0("GO-",name))
    create_dir(save_folder_GO)
    GO_analysis(gene_list,name,save_folder_GO)
    # KEGG analysis
    save_folder_KEGG <- file.path(save_dir,paste0("KEGG-",name))
    create_dir(save_folder_KEGG)
    KEGG(gene_list,name,save_folder_KEGG)
    # print run condition
    print_color_note_type4(paste0(name," KEGG-GO DONE!!!!"))
  }
}
##----------------------------------###---------------------------------------##
