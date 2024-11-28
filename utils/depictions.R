
library(ggplot2)

# Depict the datasets with negative hessian 

indeces_negative <- AvLik_NM_Pois$negative_hessian #AvLik_NM_Pois as output of compute_analysis of AvLik_POis_NM files
indeces_similar <- c(129, 132, 139)

all_incidence_data <- read.csv("../Data/incidence_data_stream_01_100runs100199.csv")

data_plot <- all_incidence_data[all_incidence_data$ID %in% c(indeces_similar, indeces_negative), ]

data_plot <- data_plot[, -c(1, 2)]
data_plot$Time <- rep(seq(1, 51, 1), 3)
data_plot$mark <- ifelse(data_plot$ID %in% indeces_negative, 'Negative hessian', 'Positive Hessian')


ggplot(data_plot, aes(x = Time, y = I )) + geom_line(aes(group = as.factor(ID), col = as.factor(ID)), linewidth = 1) + geom_point(aes(shape = as.factor(mark)), size = 2)



#depict the datasets that are discordandt between the optim and numDeriv::hessian for 5 runs in AvLik_POis_NM

#5runs 
#those negative in numDeriv: 
hessian_Pois_5run <- read.csv('../Data/AvLik_Pois_NM_hessian_rep5.csv')
indeces_nd <- hessian_Pois_5run[which(hessian_Pois_5run$Hessian < 0), 'ID']
indeces_optim <- AvLik_NM_5_runs[which(AvLik_NM_5_runs$hessian < 0), 'ID']

all_indeces <- c(indeces_nd, indeces_optim)

all_indeces <- sort(all_indeces)
#get the optimal parmaters 
optimal_parms <- AvLik_NM_5_runs[AvLik_NM_5_runs$ID %in% all_indeces, 'input.parm']


comparative_data <- read.csv('../Data/AvLik_Pois_NM_comparativehessian_rep5.csv')

comparative_data <- comparative_data[, -c(1)]
comparative_data$mark <- ifelse(comparative_data$ID %in% indeces_nd, 'NHND', 'NHO')

colors <- as.factor(comparative_data$ID)
vlines <- data.frame('parms' = optimal_parms, 'colors' = levels(colors))

ggplot(data = comparative_data, aes(x = input.parm, y = AvLik)) + geom_line(aes(group = as.factor(ID), col = as.factor(ID)), linewidth = 2) + geom_point(aes(shape = as.factor(mark)), size = 2) + geom_vline(xintercept = optimal_parms, linetype = 2, linewidth = 1.5, col = 'red') + geom_text(x=-2.29, y= 90, label="103") + geom_text(x=-2.21, y= 90, label="141") + geom_text(x=-2.01, y= 90, label="126") + geom_text(x=-1.99, y= 90, label="112")


#depict: closer look to the functions with bad hessian in 5 runs 
data_finer_grid_bad_hessian <-  read.csv('../Data//AvLik_Pois_NM_comparativehessian_finergrid_rep5.csv' )
data_103 <- data_finer_grid_bad_hessian[data_finer_grid_bad_hessian$ID == 103, ]
data_103L <- data_103[data_103$Code == 'L', ]
plot(data_103L$input.parm, data_103L$AvLik, type = 'b')
