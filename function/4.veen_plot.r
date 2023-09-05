# Function1:CREATE dir
create_dir <- function(list_dir){
  for (i in list_dir) {
    if(! dir.exists(i)){
      dir.create(i)
    }
  }
  cat("create dir finish!!!!","\n")
}
# Function2:print run condition style 2
print_color_note_type1 <- function(sologo){
  cat(underline(bold(bgCyan("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                           "),sologo,"\n")
  cat(rep(" ",15),factor)
  cat("\n")
}
# Function3:print run condition style 3
print_color_note_type2 <- function(sologo){
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                           "),sologo,"\n")
  cat(rep(" ",15),factor)
  cat("\n")
  cat(underline(bold(bgCyan("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
}
# function2: gene symbol convert
gene_symbol_convert <- function(data_up_gene_1){
  if (data_1_tax_id != data_2_tax_id){
    temp <- homologene(data_up_gene_1, inTax = data_1_tax_id, outTax = data_2_tax_id)
    data_up_gene_1 <- temp$`10090`
  }else{
    data_up_gene_1 <- data_up_gene_1
  }
  return(data_up_gene_1)
}
# function3:up data deal
data_up_deal <- function(){
  # get data 1 up gene data
  data_1_up_gene <- subset(data_1,data_1$avg_log2FC > LFC  & data_1$p_val < Pvalue)
  # get data 2 up gene data
  data_2_up_gene <- subset(data_2,data_2$avg_log2FC > LFC  & data_2$p_val < Pvalue)
  # get inter up gene
  data_up_gene_1 <- data_1_up_gene$gene
  # convert gene symbol
  data_up_gene_1 <- gene_symbol_convert(data_up_gene_1)
  # get inter up gene
  data_up_gene_2 <- data_2_up_gene$gene
  # get up list
  up_list <- list(data_1=data_up_gene_1,data_2=data_up_gene_2)
  # get inter gene list
  up_inter_gene_list <-  intersect(data_up_gene_1,data_up_gene_2)
  up_inter_gene_list <- as.data.frame(up_inter_gene_list)
  colnames(up_inter_gene_list)[1] <- "gene_symbol"
  # save data
  write.csv(up_inter_gene_list,file.path(output_dir,paste0("LFC > ",LFC,"up","-",name_list[1],"-",name_list[2],".csv")))
  names(up_list) <- name_list
  return(up_list)
}
# function4:down data deal
data_down_deal <- function(){
  # get data 1 down gene data
  data_1_down_gene <- subset(data_1,data_1$avg_log2FC < -LFC & data_1$p_val < Pvalue)
  # get data 2 down gene data
  data_2_down_gene <- subset(data_2,data_2$avg_log2FC < -LFC & data_2$p_val < Pvalue)
  # get inter up gene
  data_down_gene_1 <- data_1_down_gene$gene
  # convert gene symbol
  data_down_gene_1 <- gene_symbol_convert(data_down_gene_1)
  # get inter up gene
  data_down_gene_2 <- data_2_down_gene$gene
  # get up list
  down_list <- list(data_1=data_down_gene_1,data_2=data_down_gene_2)
  # get inter gene list
  down_inter_gene_list <-  intersect(data_down_gene_1,data_down_gene_2)
  down_inter_gene_list <- as.data.frame(down_inter_gene_list)
  colnames(down_inter_gene_list)[1] <- "gene_symbol"
  # save data
  write.csv(down_inter_gene_list,file.path(output_dir,paste0("LFC < ",-LFC," down","-",name_list[1],"-",name_list[2],".csv")))
  names(down_list) <- name_list
  return(down_list)
}
# function5: draw veen plot
draw_veen_plot <- function(data,name){
  # data <- up_list
  p <- ggVennDiagram(data,
                     set_color = c("red","blue"),
                     set_size = 6,
                     label = "both",
                     label_percent_digit = 1,
                     label_color = "black",
                     label_alpha = 0.8,
                     label_size = 6,
                     edge_color = c("red","blue")) +
    scale_fill_distiller(palette = "Set3") +
    theme(legend.position = "none")+
    annotate("text", x=500, y=180, size = 4,colour="black",label= paste0(name," Gene Cutoff:p_val < ",Pvalue ," & LFC > ",LFC))
  # save plot
  ggsave(filename = file.path(figure_dir,paste0("LFC > ",LFC,name,"-",name_list[1],"-",name_list[2],".pdf")),device = "pdf",height = 4,width = 5,plot = p)
  ggsave(filename = file.path(figure_dir,paste0("LFC > ",LFC,name,"-",name_list[1],"-",name_list[2],".png")),device = "png",height = 4,width = 5,plot = p,dpi=1200)
}
