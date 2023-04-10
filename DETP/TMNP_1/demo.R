#load function from helper_func.R and helper_plot.R

# we use RNA-seq diffierence analysis result DH_DEG.csv  as demo
DH_DEG <- read.csv("DH_DEG.csv",header = T,row.names = 1)
#predict
DH_target_result <- fun_target(DH_DEG)
DH_target_result <- DH_target_result$query_target_gsea
write.csv(DH_target_result,file = "../../result/DH_target_result.csv")

res_name <- colnames(DH_DEG)#res_name represent the group 

plot_fun(DH_target_result,Module=1,res_name)

a = DH_target_result$`Target name`
write.table(a,file = "a.txt",col.names = F,row.names = F)
