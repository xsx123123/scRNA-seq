# author:zhang jian
# date:2023.8.17
# version:1.0v
# function1:create dir
# PATH-2 loading function
# function1:create dir
create_dir <- function(list_dir){
  for (i in list_dir) {
    if(! dir.exists(i)){
      dir.create(i)
    }
  }
}
# function2:print run condition UP
print_color_note_UP <- function(logo){
  # input print with print
  if (is.na(width_print)){
    width_print <- 100
  }else{
    width_print <- width_print
  }
  # print subfactor of print dafult:(●´∀｀●)ﾉ
  logo_style_1 <- "(●´∀｀●)ﾉ"
  logo_style_2 <- "(￣▽￣)~*"
  padded_text_style1 <- str_pad(logo_style_1, width=width_print, pad = "#", side = "both")
  padded_text_style2 <- str_pad(logo_style_2, width=width_print, pad = "#", side = "both")
  # pint
  cat(bold(cyan(padded_text_style1)),"\n")
  cat("\n")
  cat(bold(str_pad(logo, width=width_print, pad = " ", side = "both"),"\n"))
  cat("\n")
}
# function3:print run condition middle
print_color_note_middle <- function(logo){
  # input print with print
  if (is.na(width_print)){
    width_print <- 100
  }else{
    width_print <- width_print
  }
  # print subfactor of print dafult:(●´∀｀●)ﾉ
  cat("\n")
  cat(bold(str_pad(logo, width=width_print, pad = " ", side = "both"),"\n"))
  cat("\n")
}
# function4:print run condition DOWN
print_color_note_DOWN <- function(logo){
  # input print with print
  if (is.na(width_print)){
    width_print <- 100
  }else{
    width_print <- width_print
  }
  # print subfactor of print dafult:(●´∀｀●)ﾉ
  logo_style_1 <- "(●´∀｀●)ﾉ"
  logo_style_2 <- "(￣▽￣)~*"
  padded_text_style1 <- str_pad(logo_style_1, width=width_print, pad = "#", side = "both")
  padded_text_style2 <- str_pad(logo_style_2, width=width_print, pad = "#", side = "both")
  # pint
  cat("\n")
  cat(bold(str_pad(logo, width=width_print, pad = " ", side = "both"),"\n"))
  cat("\n")
  cat(bold(cyan(padded_text_style2)))
  cat("\n")
}
# function5:print run condition WARRING
print_color_note_warring <- function(logo){
  # input print with print
  if (is.na(width_print)){
    width_print <- 100
  }else{
    width_print <- width_print
  }
  # print subfactor of print dafult:(●´∀｀●)ﾉ
  logo_style_1 <- "(σ｀д′)σ"
  logo_style_2 <- "(σ｀д′)σ"
  padded_text_style1 <- str_pad(logo_style_1, width=width_print, pad = "#", side = "both")
  padded_text_style2 <- str_pad(logo_style_2, width=width_print, pad = "#", side = "both")
  # pint
  cat(bold(bgRed(padded_text_style1)),"\n")
  cat("\n")
  cat(bold(str_pad(logo, width=width_print, pad = " ", side = "both"),"\n"))
  cat("\n")
  cat(bold(bgRed(padded_text_style2)))
  cat("\n")
}
# function6:print run condition NOTE
print_color_note <- function(logo){
  # input print with print
  if (is.na(width_print)){
    width_print <- 100
  }else{
    width_print <- width_print
  }
  # print subfactor of print dafult:(●´∀｀●)ﾉ
  logo_style_1 <- "(●´∀｀●)ﾉ"
  logo_style_2 <- "(￣▽￣)~*"
  padded_text_style1 <- str_pad(logo_style_1, width=width_print, pad = "#", side = "both")
  padded_text_style2 <- str_pad(logo_style_2, width=width_print, pad = "#", side = "both")
  # pint
  cat(bold(cyan(padded_text_style1)),"\n")
  cat("\n")
  cat(bold(str_pad(logo, width=width_print, pad = " ", side = "both"),"\n"))
  cat("\n")
  cat(bold(cyan(padded_text_style2)))
  cat("\n")
}
# function7:draw QC plot of raw data
sc_RNA_seq_raw_qc <- function(data,name){
  # data <- qc_data
  # print run condition
  cat("\n")
  cat(cyan("DRAW scRNA-seq RAW data plot DO!!!!!!!"))
  cat("\n")
  # draw plot
  p <- VlnPlot(data,group.by = "orig.ident", 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3,
               pt.size=0.05,
               combine = TRUE)
  # plot figure
  print(p)
  # save plot pdf png format
  ggsave(file.path(figure_dir,paste0("1.",name," raw data qc.pdf")),plot = p,width = 20,height = 20,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,paste0("1.",name," raw data qc.png")),plot = p,width = 20,height = 20,units = "cm",device = "png",dpi=1000)
  # print run condition
  cat("\n")
  cat(cyan("DRAW scRNA-seq RAW data plot DONE!!!!!!!","\n"))
  cat("\n")
}
# function8:draw QC plot of clean data
sc_RNA_seq_filted_qc <- function(data,name){
  # print run condition
  cat("\n")
  cat(cyan("DRAW scRNA-seq filted data plot DO!!!!!!!"))
  cat("\n")
  # draw plot
  p <- VlnPlot(data,group.by = "orig.ident", 
               features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
               ncol = 3,
               pt.size=0.05,
               combine = TRUE)
  print(p)
  # save plot pdf png format
  ggsave(file.path(figure_dir,paste0("2.",name," filted data qc.pdf")),plot = p,width = 20,height = 20,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,paste0("2.",name," filted data qc.png")),plot = p,width = 20,height = 20,units = "cm",device = "png",dpi=1000)
  # print run condition
  cat("\n")
  cat(cyan("DRAW scRNA-seq filted data plot DONE!!!!!!!","\n"))
  cat("\n")
}
# function9:input Seaurt & qc data
merge_data <- function(dir){
  scrna_seq_list <- list()
  for(i in 1:length(dir)){
    counts <- Read10X(data.dir = dir[i])
    scrna_seq_list[[i]] <- CreateSeuratObject(counts, min.cells=1)
    scrna_seq_list[[i]][["percent.mt"]] <- PercentageFeatureSet(scrna_seq_list[[i]], pattern = "^[AM][tY]")
    # draw raw data QC result
    sc_RNA_seq_raw_qc(scrna_seq_list[[i]],names(dir[i]))
    # filted data
    scrna_seq_list[[i]] <- subset(scrna_seq_list[[i]], subset = nFeature_RNA > nFeature_RNA_cutoff_1 & nFeature_RNA < nFeature_RNA_cutoff_2 & percent.mt < percent.mt_cutoff)
    # draw clean data QC result
    sc_RNA_seq_filted_qc(scrna_seq_list[[i]],names(dir[i]))
    # print run condition
    cat("\n")
    cat(paste0("Input ",names(dir[i])," scRNA-seq data !!!"),"\n")
  }
  return(scrna_seq_list)
}
# function10:normal & find Variable Features
normal_fing_vari_feature <- function(data){
  #data <- scrna_seq_list
  for (i in 1:length(data)) {
    data[[i]] <- NormalizeData(data[[i]],verbose = F)
    data[[i]] <- FindVariableFeatures(data[[i]], selection.method = "vst",nfeatures = 2000, verbose = F)
  }
  return(data)
}
# function11:draw UMAP plot
dimplot_double <- function(data,name,figure_num){
  # data <- scrna_seq_integr
  p <- DimPlot(data, reduction = "umap", split.by = "orig.ident",label = TRUE,alpha=1,pt.size=0.4,label.size=3.2)+
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(name," UMAP plot")) +
    labs(x="UMAP-1",y="UMAP-2")+
    scale_y_continuous(n.breaks = 20) +
    scale_x_continuous(n.breaks = 20) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="bottom",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),    
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))  +
    guides(color=guide_legend(nrow=2,title="Cluster-ID",title.position = "left",override.aes = list(size=3,alpha=1)))
  # save
  ggsave(file.path(figure_dir,paste0(figure_num," ",name," scRNA-seq-UMAP-spilt-plot.pdf")),plot=p,width = 20,height = 14,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,paste0(figure_num," ",name," scRNA-seq-UMAP-spilt-plot.png")),plot=p,width = 20,height = 14,units = "cm",device = "png",dpi=1000)
}
# function12:draw UMAP plot of integr data
integr_DimPlot <- function(data,name){
  # name <- "Integr"
  # data <- scrna_seq_integr
  if (length(levels(data@active.ident)) < 7){
    nrow = 1
    hei = 13
  }else{
    if (length(levels(data@active.ident)) < 14){
      nrow = 2
      hei = 14
    }else{
      if (length(levels(data@active.ident)) < 21){
        nrow = 3
        hei = 15
      }else{
        nrow = 4
        hei = 16
      }
    }}
  # draw integr DimPlot
  p <- DimPlot(data, reduction = "umap",group.by='orig.ident')+
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(name," UMAP plot")) +
    labs(x="UMAP-1",y="UMAP-2")+
    scale_y_continuous(n.breaks = 20) +
    scale_x_continuous(n.breaks = 20) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="bottom",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),    
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))+
    guides(color=guide_legend(nrow=nrow,title="Cluster-ID",title.position = "left",override.aes = list(size=2.5,alpha=0.8)))
  # ggsave
  ggsave(file.path(figure_dir,"5.integr_group_UMAP-plot.pdf"),plot=p,width = 15,height = 16,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,"5.integr_group_UMAP-plot.png"),plot=p,width = 15,height = 16,units = "cm",device = "png",dpi=1000)
  # draw integr DimPlot cell cluster
  p1 <- DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE) +
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(name," UMAP plot")) +
    labs(x="UMAP-1",y="UMAP-2")+
    scale_y_continuous(n.breaks = 20) +
    scale_x_continuous(n.breaks = 20) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="bottom",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),    
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))+
    guides(color=guide_legend(nrow=nrow,title="Cluster-ID",title.position = "left",override.aes = list(size=2.5,alpha=0.8)))
  ggsave(file.path(figure_dir,"6.integr_group_UMAP-plot.pdf"),plot=p1,width = 15,height = 16,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,"6.integr_group_UMAP-plot.png"),plot=p1,width = 15,height = 16,units = "cm",device = "png",dpi=1000)
  # merge two plot
  all_p <- ggarrange(p,p1,ncol=2,nrow=1,labels=c("A","B"))
  ggsave(file.path(figure_dir,"7.integr_scRNA-seq-UMAP-plot.pdf"),plot=all_p,width = 26,height = hei,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,"7.integr_scRNA-seq-UMAP-plot.png"),plot=all_p,width = 26,height = hei,units = "cm",device = "png",dpi=1000)
}
# function13:draw UMAP plot of integr of two data
dimplot_double_2 <- function(data,name){
  # data <- scrna_seq_integr
  p <- DimPlot(data, reduction = "umap", split.by = "orig.ident",label = TRUE,alpha=1,pt.size=0.4,label.size=3.2)+
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(name," UMAP plot")) +
    labs(x="UMAP-1",y="UMAP-2")+
    scale_y_continuous(n.breaks = 20) +
    scale_x_continuous(n.breaks = 20) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="bottom",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),    
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))  +
    guides(color=guide_legend(nrow=2,title="Cluster-ID",title.position = "left",override.aes = list(size=3,alpha=1)))
  # save
  ggsave(file.path(figure_dir,"8.integr_scRNA-seq-UMAP-spilt-plot.pdf"),plot=p,width = 20,height = 14,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,"8.integr_scRNA-seq-UMAP-spilt-plot.png"),plot=p,width = 20,height = 14,units = "cm",device = "png",dpi=1000)
}
# function14:find marker gene for everyone cell cluster(normal)
get_cell_cluster_markder_gene <- function(data){
  if(threads == T){
    cat("NOT USE THIS PROCESS (PS:SPEED VERY SLOW)","\n")
  }else{cluster_number <- levels(data@meta.data$seurat_clusters)
  all_cluster_top10_gene <- data.frame()
  # process 
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total = length(cluster_number), clear = FALSE, width = 100
  )
  # get everyone cell cluster marker gene
  for(i in cluster_number){
    temp_cell_cluster_marker_gene <- FindConservedMarkers(data, ident.1 = i, grouping.var = "orig.ident", verbose = FALSE)
    temp_cell_cluster_marker_gene$cell_cluster <- paste0("cell cluster ",i)
    # save cell cluster all marker gene
    write.csv(temp_cell_cluster_marker_gene,file.path(marker_gene_output_dir,paste0("Cell cluster ",i,"marker gene list.csv")))
    # get top10 cell cluster
    temp_cell_cluster_marker_gene_top10 <- head(temp_cell_cluster_marker_gene,10)
    all_cluster_top10_gene <- rbind(all_cluster_top10_gene,temp_cell_cluster_marker_gene_top10)
    pb$tick()
  }
  # save
  write.csv( all_cluster_top10_gene,file.path(marker_gene_output_dir,paste0("All marker top 10 gene list.csv")))}
  
}
# function15:find marker gene for everyone cell cluster(multithreading)
get_cell_cluster_markder_gene_threads <- function(i){
  # get everyone cell cluster marker gene
  temp_cell_cluster_marker_gene <- FindConservedMarkers(scrna_seq_integr, ident.1 = i, grouping.var = "orig.ident", verbose = FALSE)
  temp_cell_cluster_marker_gene$cell_cluster <- paste0("cell cluster ",i)
  # save cell cluster all marker gene
  write.csv(temp_cell_cluster_marker_gene,file.path(marker_gene_output_dir,paste0("Cell cluster ",i," marker gene list.csv")))
  # print run condition
}
# function16:(multithreading)
multithreading_get_marker_gene <- function(){
  if (threads == T){
    cluster_number <- levels(scrna_seq_integr@meta.data$seurat_clusters)
    # use threads number
    num_cores = detectCores() - 4
    # print run condition
    print_color_note_UP("Use multithreading for find marker gene DO!!!!!")
    print_color_note_middle(paste0("Use ",num_cores," Core!!!!"))
    time <- now()
    print_color_note_middle(time)
    # call multithreading
    mclapply(cluster_number, get_cell_cluster_markder_gene_threads, mc.cores = num_cores)
    # get TOP10 marker gene for everyone cell cluster
    time <- now()
    print_color_note_middle(time)
    print_color_note_DOWN("Use multithreading for find marker gene DONE!!!!!")
  }else{
    print_color_note("Don't use multithreading find marker gene")
  }
  
}
# function17:extert TOP 10 marker gene for cell cluster
Extert_top10_marker_gene <- function(list){
  top10_marker_gene <- data.frame()
  for (i in list){
    temp <- read.csv(file.path(marker_gene_output_dir,i))
    temp_top10_marker_gene <- head(temp,10)
    top10_marker_gene <- rbind(top10_marker_gene,temp_top10_marker_gene)
  }
  write.csv(top10_marker_gene,file.path(marker_gene_output_dir,"top10_marker_gene.csv"),row.names=F)
}
# function18:draw UMAP  plot
draw_umap_plot <- function(data){
  # print run condition
  cat("\n")
  cat(cyan("PCA dim PLOT (UMAP) DO!!!!!"))
  cat("\n")
  p <- DimPlot(data, reduction = "umap") 
  p1 <- p +
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(project_name," UMAP plot")) +
    labs(x="UMAP-1",y="UMAP-2")+
    scale_y_continuous(n.breaks = 20) +
    scale_x_continuous(n.breaks = 20) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="bottom",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),    
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))+
    guides(color=guide_legend(nrow=2,title="Cluster-ID",title.position = "left",override.aes = list(size=2.5,alpha=0.8)))
  # ggsave
  # print run condition
  cat("\n")
  cat(cyan("PCA dim PLOT (UMAP) DONE!!!!!"))
  cat("\n")
  return(p1)
}
# function19:draw tsne  plot
draw_tsne_plot <- function(data){
  # print run condition
  cat("\n")
  cat(cyan("PCA dim PLOT (tSNE) DO!!!!!"))
  cat("\n")
  p <- DimPlot(data, reduction = "tsne") 
  p1 <- p +
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(project_name," tSNE plot")) +
    labs(x="tSNE-1",y="tSNE-2")+
    scale_y_continuous(n.breaks = 20) +
    scale_x_continuous(n.breaks = 20) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="bottom",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),    
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))+
    guides(color=guide_legend(nrow=2,title="Cluster-ID",title.position = "left",override.aes = list(size=2.5,alpha=0.8)))
  # ggsave
  # print run condition
  cat("\n")
  cat(cyan("PCA dim PLOT (tSNE) DONE!!!!!"))
  cat("\n")
  return(p1)
}
# function20:manual annotation cell cluster
manual_annotation_figure <- function(data,genes_to_check){
  # print run condition
  cat("\n")
  cat(cyan(paste0("manual annotation ",genes_to_check," PLOT  DO!!!!!")))
  cat("\n")
  # draw  DotPlot
  p1 <- DotPlot(data, features = genes_to_check,assay='RNA')+
    theme_classic()+
    ggtitle(paste0(project_name," DotPlot for  special gene")) +
    labs(x="Feature",y="cell cluster ID")+
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 12),
          legend.position="right",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10,face="bold"),
          axis.title.y = element_text(size = 10,face="bold"),    
          axis.text.x = element_text(size = 8,face="bold"),
          axis.text.y = element_text(size = 8,face="bold"),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))
  # draw VlnPlot
  p2 <- VlnPlot(data, features = genes_to_check)+
    theme_classic()+
    labs(x="cell cluster ID",y="Expression level")+
    scale_y_continuous(expand = c(0,0),n.breaks = 6) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="none",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10,face="bold"),
          axis.title.y = element_text(size = 10,face="bold"),    
          axis.text.x = element_text(size = 8,face="bold"),
          axis.text.y = element_text(size = 8,face="bold"),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))
  # draw FeaturePlot
  p3 <- FeaturePlot(data, features = genes_to_check)+
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(project_name," ",genes_to_check," UMAP plot")) +
    labs(x="UMAP-1",y="UMAP-2")+
    scale_y_continuous(n.breaks = 20) +
    scale_x_continuous(n.breaks = 20) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 15),
          legend.position="bottom",
          legend.title=element_text(size = 10), 
          text = element_text(size = 8,family="sans"),           
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),    
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))+
    guides(color=guide_legend(nrow=2,title="Cluster-ID",title.position = "left",override.aes = list(size=2.5,alpha=0.8)))
  # merge plot
  all_p <- ggarrange(umap_plot,p3,p1,p2,ncol=2,nrow=2,labels=c("A","B","C","D"))
  # print plot
  print(all_p)
  # save plot
  ggsave(file.path(annotation_figure_dir,paste0(genes_to_check," merge_plot.pdf")),device = "pdf",width = 10,height = 12,plot=all_p)
  ggsave(file.path(annotation_figure_dir,paste0(genes_to_check," merge_plot.png")),device = "png",width = 10,height = 12,plot=all_p,dpi=600)
  # print run condition
  cat("\n")
  cat(cyan(paste0("manual annotation ",genes_to_check," PLOT  DONE!!!!!")))
  cat("\n")
}
# function21: annotation_cell_cluster_name tsne
annotation_cell_cluster_name <- function(data,new.cluster.ids){
  #data <- scrna_seq_integr
  names(new.cluster.ids) <- levels(data)
  rename_data <- RenameIdents(data, new.cluster.ids)
  p <-  DimPlot(rename_data,split.by = "orig.ident", reduction = "tsne",label = TRUE, pt.size = 0.25,label.size=1.5) + NoLegend() +
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(project_name," tSNE plot")) +
    labs(x="tSNE-1",y="tSNE-2")+
    scale_y_continuous(n.breaks = 10) +
    scale_x_continuous(n.breaks = 10) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 8),
          legend.position="right",
          legend.title=element_text(size = 6,face="bold"), 
          text = element_text(size = 4,face="bold"),           
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),    
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))+
    guides(color=guide_legend(title="Cluster-ID",title.position = "top",override.aes = list(size=2,alpha=0.8)))
  print(p)
  # save
  ggsave(file.path(figure_dir,paste0("9.tSNE-annotation-plot.pdf")),plot=p,width = 18,height = 8,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,paste0("9.tSNE-annotation-plot.png")),plot=p,width = 18,height = 8,units = "cm",device = "png",dpi=1000)
  # print run condition
  p <-  DimPlot(rename_data, split.by = "orig.ident",reduction = "umap", label = TRUE, pt.size = 0.25,label.size=1.5) + NoLegend() +
    geom_vline(xintercept=0,lty=4,col="black",lwd=0.3) +  
    geom_hline(yintercept = 0,lty=4,col="black",lwd=0.3) +
    theme_classic()+
    ggtitle(paste0(project_name," UMAP plot")) +
    labs(x="UMAP-1",y="UMAP-2")+
    scale_y_continuous(n.breaks = 10) +
    scale_x_continuous(n.breaks = 10) +
    theme(panel.border = element_rect(linetype = 1, fill = NA,size=0.6),
          plot.title = element_text(hjust = 0.5,size = 8),
          legend.position="right",
          legend.title=element_text(size = 6,face="bold"), 
          text = element_text(size = 4,face="bold"),           
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),    
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))+
    guides(color=guide_legend(title="Cluster-ID",title.position = "top",override.aes = list(size=2,alpha=0.8)))
  print(p)
  # save
  ggsave(file.path(figure_dir,paste0("9.UMAP-annotation-plot.pdf")),plot=p,width = 18,height = 8,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,paste0("9.UMAP-annotation-plot.png")),plot=p,width = 18,height = 8,units = "cm",device = "png",dpi=1000)
  # print run condition
  return(rename_data)
}
# function22: annotation_cell_cluster_name doplot
draw_annotation_cell_cluster_doplot <- function(data,genes_to_check){
  #data <- rename_maritx
  p<- DotPlot(data, features = genes_to_check,assay='RNA',dot.scale=2.5) + 
    coord_flip()+
    theme_classic()+
    ggtitle(paste0(project_name," Dimplot")) +
    labs(x="",y="")+
    theme(plot.title = element_text(hjust = 0.5,size = 8),
          legend.position="right",
          legend.title=element_text(size = 5), 
          text = element_text(size = 9),title = element_text(size = 7),           
          axis.title.x = element_text(size = 8),axis.title.y = element_text(size = 7),    
          axis.text.x = element_text(size = 5,angle=45,hjust=0.9,vjust=0.9),axis.text.y = element_text(size = 5),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="black",size=0.2,lineend = 1),
          axis.ticks.y=element_line(color="black",size=0.2,lineend = 1))
  print(p)
  # save plot
  # save
  ggsave(file.path(figure_dir,paste0("10.annotation-DotPlot.pdf")),plot=p,width = 8,height = 11,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,paste0("10.annotation-DotPlot.png")),plot=p,width = 8,height = 11,units = "cm",device = "png",dpi=1000)
  # print run condition
}
# function23: states cell cluster rate
states_cell_cluster_rate <- function(data){
  #data <-rename_maritx
  if (length(levels(data@active.ident)) <= 10){
    ncol = 1
    width = 8
  }else{
    if (length(levels(data@active.ident)) <= 20){
      ncol = 2
      width = 11
    }else{
      if (length(levels(data@active.ident)) <= 30){
        ncol = 3
        width = 14
      }
    }
  }
  # extert cell cluster number
  rate_data <- as.data.frame(table(Idents(data),data$orig.ident))
  # convert rate
  rate_data_1 <- data.frame()
  for(i in levels(rate_data$Var2)){
    temp_data <- subset(rate_data,rate_data$Var2==i)
    temp_sum <- sum(temp_data$Freq)
    temp_data$rate <- temp_data$Freq/temp_sum
    rate_data_1 <- rbind(rate_data_1,temp_data)
  }
  # convert %
  rate_data_1$rate <- rate_data_1$rate*100
  # change col name
  colnames(rate_data_1)[1] <- "Cell_cluster_name"
  colnames(rate_data_1)[2] <- "group"
  colnames(rate_data_1)[3] <- "frequently"
  colnames(rate_data_1)[4] <- "rate"
  # draw plot
  colour1 <- c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  p <- ggplot(rate_data_1,aes(x=group,y=rate,fill=Cell_cluster_name))+
    geom_bar(position = 'stack',stat="identity",width=0.4)+
    labs(x="group name",y = "Cell cluster frequently (%)",title = "WT Cell cluster rate")+
    scale_fill_manual(values=colour1)+
    scale_y_continuous(expand = c(0,0),n.breaks = 10) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size =6,family="sans",face = "bold"),
          legend.title=element_text(size =6,face = "bold",family="sans"), 
          legend.text =element_text(size =4,face = "bold",family="sans"), 
          text = element_text(size = 4,face = "bold"),title = element_text(size = 4),           
          axis.title.x = element_text(size = 5,family="sans",face = "bold"),axis.title.y = element_text(size = 5,family="sans",face = "bold"),    
          axis.text.x = element_text(size = 4,family="sans"),axis.text.y = element_text(size = 4,family="sans"),
          panel.border = element_rect(color = "#606c70", fill = NA, size = 0.3),
          axis.line.x=element_line(linetype=1,color="#606c70",size=0.12),       
          axis.line.y=element_line(linetype=1,color="#606c70",size=0.12),
          axis.ticks.x=element_line(color="#606c70",size=0.12,lineend = 0.05),
          axis.ticks.length=unit(.08,"lines"),
          axis.ticks.y=element_line(color="#606c70",size=0.12,lineend = 0.05))+
    guides(fill=guide_legend(ncol=ncol,title="Cell cluster ID",override.aes = list(size=0.02,alpha=1)))
  print(p)
  # save
  write.csv(rate_data_1,file.path(output_dir,paste0("2.proportion_annotation-plot.csv")),row.names = F)
  ggsave(file.path(figure_dir,paste0("11.proportion_annotation-plot.pdf")),plot=p,width = width,height = 8,units = "cm",device = "pdf")
  ggsave(file.path(figure_dir,paste0("11.proportion_annotation-plot.png")),plot=p,width = width,height = 8,units = "cm",device = "png",dpi=1000)
  # print run condition
}
# function24: DRAW volcano PLOT
draw_volcano <- function(deg,i){
  # draw plot
  deg_result <- deg
  deg_result$log10 <- -log10(deg_result$p_val)
  # ADD UP&DOWN&NO GENE TAG
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$p_val < 0.05) & (deg_result$avg_log2FC > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$p_val  < 0.05) & (deg_result$avg_log2FC < -1))] = "Down-regulated"
  # SORT BY p_val VALUE
  deg_result <- deg_result[order(deg_result$p_val),]
  # GET NOT&UP&DOWN DATA
  non_deg_result <- subset(deg_result,deg_result$Group =="Non-significan")
  up_deg_result <- subset(deg_result,deg_result$Group =="Up-regulated")
  down_deg_result <- subset(deg_result,deg_result$Group =="Down-regulated")
  # GET TOP 15 p_val GENE UP&DOWN
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),15)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),15)
  # get y_aes
  y_aes <- deg_result$log10
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$avg_log2FC
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  p <- ggplot(deg_result, aes(x = avg_log2FC, y = log10)) + 
    geom_point(data=non_deg_result,aes(x = avg_log2FC, y = log10),size=0.02,shape = 21,color="#C7C7C7",alpha=0.25) +
    geom_point(data=deg_result_up,aes(x = avg_log2FC, y = log10),size=0.02,shape = 21,fill="#e41749",alpha=0.5) +
    geom_point(data=up_deg_result,aes(x = avg_log2FC, y = log10),size=0.02,shape = 21,color="#e41749",alpha=0.4) +
    geom_point(data=deg_result_down,aes(x = avg_log2FC, y = log10),size=0.02,shape = 21,fill="#41b6e6",alpha=0.5) +
    geom_point(data=down_deg_result,aes(x = avg_log2FC, y = log10),size=0.02,shape = 21,color="#41b6e6",alpha=0.4) +
    geom_vline(xintercept=0,lty=2,col="black",lwd=0.1) +  
    geom_hline(yintercept = 1.3,lty=2,col="black",lwd=0.1) +  
    labs(x= bquote("RNA-seq " * log[2] * " fold change " * .(EXP_NAEE) * ""),y= expression(paste(-log[10], "P-value")),title =paste0(i," Volcano Plot")) + 
    geom_text_repel(data = deg_result_up,aes(avg_log2FC, log10, label= gene),size=1,colour="black",fontface="bold.italic",
                    segment.alpha = 0.5,segment.size = 0.15,segment.color = "black",min.segment.length=0,
                    box.padding=unit(0.2, "lines"),point.padding=unit(0, "lines"),force = 20,max.iter = 3e3,
                    xlim=c(0, 6),max.overlaps = 25,arrow=arrow(length = unit(0.02, "inches"))) + 
    geom_text_repel(data = deg_result_down,aes(avg_log2FC, log10, label= gene),size=1,colour="black",fontface="bold.italic",
                    segment.alpha =0.5,segment.size = 0.15,segment.color = "black",min.segment.length=0,
                    box.padding=unit(0.2, "lines"),point.padding=unit(0, "lines"),force = 20,max.iter = 3e3,
                    xlim=c(-6, 0),max.overlaps = 25,arrow=arrow(length = unit(0.02, "inches"))) + 
    annotate("text", x=-(x_aes*0.8), y=y_aes_value*0.9,size = 1.5,colour="#41b6e6", label= paste0("Down in ",EXP_NAEE))+
    geom_segment(aes(x = -(x_aes*0.5), y = y_aes_value*0.85, xend = -(x_aes*1.05), yend = y_aes_value*0.85),arrow = arrow(length = unit(0.15, "cm")),colour="#41b6e6",linewidth=0.4,
                 alpha=0.15,lineend="round",linejoin="round") +
    annotate("text", x=x_aes*0.8, y=y_aes_value*0.9, size = 1.5,colour="#e41749",label= paste0("UP in ",EXP_NAEE) )+
    geom_segment(aes(x = (x_aes*0.5), y = y_aes_value*0.85, xend = (x_aes*1.05), yend = y_aes_value*0.85),arrow = arrow(length = unit(0.15, "cm")),colour="#e41749",linewidth=0.4,
                 alpha=0.15,lineend="round",linejoin="round") +
    guides(color=guide_legend(override.aes = list(size=10)),) + 
    scale_x_continuous(limits=c(-(x_aes*1.2),(x_aes*1.2)),n.breaks = 8) +
    scale_y_continuous(limits=c(0,y_aes_value),n.breaks = 10) +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5,size =5,family="sans",face = "bold"),legend.position="none",
          legend.title=element_text(size =4,face = "bold",family="sans"), 
          text = element_text(size = 4,family="sans"),title = element_text(size = 4),           
          axis.title.x = element_text(size = 4,family="sans",face = "bold"),axis.title.y = element_text(size = 4,family="sans",face = "bold"),    
          axis.text.x = element_text(size = 4,family="sans"),axis.text.y = element_text(size = 4,family="sans"),
          panel.border = element_rect(color = "#606c70", fill = NA, size = 0.3),
          axis.line.x=element_line(linetype=1,color="#606c70",size=0.12),       
          axis.line.y=element_line(linetype=1,color="#606c70",size=0.12),
          axis.ticks.x=element_line(color="#606c70",size=0.12,lineend = 0.05),
          axis.ticks.length=unit(.08,"lines"),
          axis.ticks.y=element_line(color="#606c70",size=0.12,lineend = 0.05))
  # save plot
  ggsave(file.path(deg_figure_dir,paste0(i," Volcano Plot.pdf")),plot = p,width = 7,height = 5,units = "cm")
  ggsave(file.path(deg_figure_dir,paste0(i," Volcano Plot.png")),device = "png",plot = p,width = 7,height = 5,units = "cm",dpi = 1000)
}
# function25: draw activate  volcano plot
draw_volcano_activate <- function(deg,i){
  # draw plot
  deg_result <- deg
  deg_result$log10 <- -log10(deg_result$p_val)
  # ADD UP&DOWN&NO GENE TAG
  deg_result$label = NA
  deg_result$Group <- "Non-significan"
  deg_result$Group[which((deg_result$p_val < 0.05) & (deg_result$avg_log2FC > 1))] = "Up-regulated"
  deg_result$Group[which((deg_result$p_val  < 0.05) & (deg_result$avg_log2FC < -1))] = "Down-regulated"
  # SORT BY p_val VALUE
  deg_result <- deg_result[order(deg_result$p_val),]
  # GET NOT&UP&DOWN DATA
  non_deg_result <- subset(deg_result,deg_result$Group =="Non-significan")
  up_deg_result <- subset(deg_result,deg_result$Group =="Up-regulated")
  down_deg_result <- subset(deg_result,deg_result$Group =="Down-regulated")
  # GET TOP 15 p_val GENE UP&DOWN
  deg_result_up <- head(subset(deg_result,deg_result$Group == "Up-regulated"),15)
  deg_result_down <- head(subset(deg_result,deg_result$Group == "Down-regulated"),15)
  # get y_aes
  y_aes <- deg_result$log10
  # get y_aes
  y_aes <- deg_result$log10
  # remove inf
  y_aes <- y_aes[is.finite(y_aes)]
  y_1 <- sort(y_aes, decreasing = TRUE)[1]
  y_2 <- sort(y_aes, decreasing = TRUE)[2]
  if ( max(y_aes) > 300){
    y_aes_value <- 250
  }else{
    if(y_1/y_2 > 1.4){
      y_aes_value <- (y_1+y_2)/2
    }else{
      y_aes_value <- max(y_aes)*1.05
    }
  }
  # x_aes deal
  x_aes <- deg_result$avg_log2FC
  # remove NA
  x_aes <- na.omit(x_aes)
  if (max(x_aes) > 7.5){
    x_aes <- 7.5
  }else{
    x_aes <- max(x_aes)
  }
  # draw volcano plot
  plot <- plot_ly(deg_result,
                  x = ~avg_log2FC,
                  y = ~-log10(p_val),
                  text = ~gene, 
                  type = 'scatter',
                  mode = 'markers',
                  alpha = 0.6,
                  marker = list(
                    color = ~I(ifelse(avg_log2FC > 0.5 & -log10(p_val) > 1.3, "red",
                                      ifelse(avg_log2FC < -0.5 & -log10(p_val) > 1.3, "blue", "gray")))))
  
  layout(plot, xaxis = list(range = c(-(max(deg_result$avg_log2FC)*1.3), 1.3 * max(deg_result$avg_log2FC))),
         yaxis = list(range = c(0, max(-log10(deg_result$p_val)))),
         title = paste0(i," volcano plot"))
  #
  saveWidget(plot, file = file.path(deg_figure_dir,paste0(i," Volcano plot.html")), selfcontained = TRUE)
}
# function26: single cell cluster DEG
DEG_single_cell_cluster <- function(i){
  print_color_note_UP(paste0("Cell cluster ",i," DEG DO!!!!"))
  rename_maritx@meta.data$celltype <- rename_maritx@active.ident
  # extert cell cluster for group1
  cells1 <- subset(rename_maritx@meta.data, celltype %in% c(i)& orig.ident == name_dir[1])  %>% rownames()
  # extert cell cluster for group2
  cells2 <- subset(rename_maritx@meta.data, celltype %in% c(i)& orig.ident == name_dir[2])  %>% rownames()
  # use MAST DEG
  deg <- FindMarkers(rename_maritx, ident.1 = cells1, ident.2 = cells2,test.use="MAST")
  deg <- data.frame(gene = rownames(deg), deg)
  # draw volcano
  draw_volcano(deg,i)
  draw_volcano_activate(deg,i)
  # set save dir name
  save_path <- file.path(deg_dir,i)
  create_dir(c(save_path))
  up_gene <- subset(deg,deg$avg_log2FC > 0.5 & deg$p_val_adj < 0.05)
  down_gene <- subset(deg,deg$avg_log2FC < -0.5 & deg$p_val_adj < 0.05)
  # save deg result
  write.csv(deg,file.path(save_path,paste0(i," DEG-result.csv")))
  write.csv(up_gene,file.path(save_path,paste0(i," UP-DEG-result.csv")))
  write.csv(down_gene,file.path(save_path,paste0(i," DOWN-DEG-result.csv")))
  # PRINT RUN CONDITION
  print_color_note_UP(paste0("Cell cluster ",i," DEG DONE!!!!"))
}
# function27: multithreading deg
multithreading_DEG <- function(){
  if (threads == T){
    cell_cluster_list <- levels(rename_maritx@active.ident)
    # use threads number
    num_cores = detectCores() - 20
    # print run condition
    print_color_note_UP("Use multithreading for DEG DO!!!!!")
    print_color_note_middle(paste0("Use ",num_cores," Core!!!!"))
    time <- now()
    print_color_note_middle(time)
    # call multithreading
    mclapply(cell_cluster_list, DEG_single_cell_cluster, mc.cores = num_cores)
    # get TOP10 marker gene for everyone cell cluster
    time <- now()
    print_color_note_middle(time)
    print_color_note_DOWN("Use multithreading for DEG DONE!!!!!")
  }else{
    print_color_note("Don't use multithreading DEG")
  }
}
# function28:print ProgressBar
ProgressBar <- function(){
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total = 10, clear = FALSE, width = 100
  )
  
  for (i in 1:10) {
    pb$tick()
    Sys.sleep(1)
  }
}
# function29:print run parameter
get_run_parameter <- function(){
  run_parameter <- c("scRNA-seq analysis root_dir is :",
                     "scRNA-seq analysis save_dir is :",
                     "scRNA-seq analysis project_name is :",
                     "scRNA-seq analysis raw data cutoff-1 is :",
                     "scRNA-seq analysis raw data cutoff-2 is :",
                     "scRNA-seq analysis raw data cutoff-3 is :",
                     "wheather use multithreading is :")
  run_parameter <- as.data.frame(run_parameter)
  run_parameter$parameter <- c(root_dir,
                               save_dir,
                               project_name,
                               nFeature_RNA_cutoff_1,
                               nFeature_RNA_cutoff_2,
                               percent.mt_cutoff,
                               threads)
  print_color_note_UP("place check run parameter!!!")
  print(run_parameter)
  print_color_note_DOWN("place check run parameter!!!")
  cat("\n")
  ProgressBar()
  cat("\n")
}
