#Case study: H2O potential energy (with PCA).
library(invgamma)
library(lhs)
library(abind)
library(mlegp)
library(laGP)
library(ggplot2)
library(coda)
library(BayesianTools)
library(readxl)
library(dplyr)
library(DiceKriging)

rm(list=ls())
x=read.csv("X_soap.csv", header=TRUE)
y=read.csv("y_energy.csv", header=TRUE)

task_id=as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(task_id+500)
source('header.R')

idx <- sample.int(nrow(x), 100)
x_train <- as.matrix(x[idx, , drop = FALSE])
x_test  <- as.matrix(x[-idx, , drop = FALSE])
y_train <- as.matrix(y[idx, , drop = FALSE])
y_test  <- as.matrix(y[-idx, , drop = FALSE])
n=100
m=21

x_train_orig <- x_train
rng_train <- apply(x_train, 2, function(col) diff(range(col)))
good_range  <- is.finite(rng_train) & (rng_train > 0)
good_finite <- colSums(!is.finite(x_train)) == 0
keep_cols   <- which(good_range & good_finite)
x_train <- x_train[, keep_cols, drop = FALSE]
x_test  <- x_test[,  keep_cols, drop = FALSE]

cat("the number of throw:", ncol(x_train_orig) - ncol(x_train),
    ", The number of keep", ncol(x_train), "\n")
#p=ncol(x_train)

max_x_train <- apply(x_train, 2, max)
min_x_train <- apply(x_train, 2, min)
scale_rng   <- pmax(max_x_train - min_x_train, .Machine$double.eps)  

xx_train <- scale(x_train, center = min_x_train, scale = scale_rng)
xx_test  <- scale(x_test,  center = min_x_train, scale = scale_rng)

pca <- prcomp(xx_train, center = FALSE, scale. = FALSE)

var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
cum_explained <- cumsum(var_explained)

k <- 50

# k <- which(cum_explained >= 0.99)[1]

Z_train <- xx_train %*% pca$rotation[, 1:k, drop = FALSE]
Z_test  <- xx_test  %*% pca$rotation[, 1:k, drop = FALSE]

prior_var=array(0,7)
prior_var[1]=1
prior_var[2]=1
prior_var[3]=1
prior_var[4]=1
prior_var[5]=4
prior_var[6]=0.5
prior_var[7]=0.5
n_out=2000
n_in=10
r=array(0,2)
r[1]=20
r[2]=5
l=0
rou=0.5

t1=Sys.time()
mlegp_model=mlegp(Z_train,y_train,constantMean = 1)
t2=Sys.time()
t_mlegp=difftime(t2,t1,units='mins')
y_pred=predict(mlegp_model,Z_test)
rmse_mlegp=sqrt(mean((y_pred-y_test)^2))/sd(y_test)

eps=sqrt(.Machine$double.eps)
t1=Sys.time()
gpi=newGPsep(Z_train,y_train,d=0.1,g=0.1*var(y_train),dK=TRUE)
mle=mleGPsep(gpi,param="both",tmin=c(eps,eps),tmax=c(10,var(y_train)))
t2=Sys.time()
t_lagp=difftime(t2,t1,units='mins')
w=predGPsep(gpi,Z_test)
rmse_lagp=sqrt(mean((w$mean-y_test)^2))/sd(y_test)

t1=Sys.time()
kriging_constant=km(~1,design=Z_train,response=y_train,covtype="gauss",nugget.estim = TRUE)
t2=Sys.time()
t_kriging_constant=difftime(t2,t1,units='mins')
y_pred=predict(kriging_constant,newdata=Z_test,type='UK',checkNames=FALSE)
rmse_kriging_constant=sqrt(mean((y_pred$mean-y_test)^2))/sd(y_test)


G=matrix(1,nrow=n,ncol=1)
#G=cbind(G,Z_train)

newG=matrix(1,nrow=m,ncol=1)
#newG=cbind(newG,Z_test)

t1=Sys.time()
sphmc_0.8=bayesian_bridge_sphmc_adp_convergence2(Z_train,y_train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=0.8,r)
t2=Sys.time()
t_sphmc_0.8=difftime(t2,t1,units='mins')

sphmc_0.8_pred=bayesian_bridge_predict(Z_train,y_train,G,Z_test,newG,sphmc_0.8)
rmse_sphmc_0.8 <- sqrt(mean((sphmc_0.8_pred-y_test)^2))/sd(y_test)

t1=Sys.time()
sphmc_1=bayesian_bridge_sphmc_adp_convergence2(Z_train,y_train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=1,r)
t2=Sys.time()
t_sphmc_1=difftime(t2,t1,units='mins')

sphmc_1_pred=bayesian_bridge_predict(Z_train,y_train,G,Z_test,newG,sphmc_1)
rmse_sphmc_1 <- sqrt(mean((sphmc_1_pred-y_test)^2))/sd(y_test)

t1=Sys.time()
sphmc_2=bayesian_bridge_sphmc_adp_convergence2(Z_train,y_train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=1.8,r)
t2=Sys.time()
t_sphmc_2=difftime(t2,t1,units='mins')

sphmc_2_pred=bayesian_bridge_predict(Z_train,y_train,G,Z_test,newG,sphmc_2)
rmse_sphmc_2 <- sqrt(mean((sphmc_2_pred-y_test)^2))/sd(y_test)

t1=Sys.time()
hmc=bayesian_bridge_sphmc_adp_hmc_convergence2(Z_train,y_train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,l=l,rou=rou)
t2=Sys.time()
t_hmc=difftime(t2,t1,units='mins')

hmc_pred=bayesian_bridge_predict(Z_train,y_train,G,Z_test,newG,hmc)
rmse_hmc <- sqrt(mean((hmc_pred-y_test)^2))/sd(y_test)


result=data.frame(task_id=task_id, rmse_mlegp=rmse_mlegp,t_mlegp=t_mlegp,rmse_lagp=rmse_lagp,t_lagp=t_lagp, rmse_kriging_constant=rmse_kriging_constant, t_kriging_constant=t_kriging_constant,rmse_sphmc_0.8=rmse_sphmc_0.8,t_sphmc_0.8=t_sphmc_0.8,rmse_sphmc_1=rmse_sphmc_1,t_sphmc_1=t_sphmc_1,rmse_sphmc_2=rmse_sphmc_2,t_sphmc_2=t_sphmc_2,rmse_hmc=rmse_hmc,t_hmc=t_hmc)

name1='Intercept'

name2=NULL
for (i in 1:k){
  name2=c(name2,bquote(Z[.(i)]))
}

new_columns_beta=data.frame(matrix(sphmc_0.8$beta,nrow=1,ncol=length(sphmc_0.8$beta)))
colnames(new_columns_beta)=name1
new_columns_omega=data.frame(matrix(sphmc_0.8$omega,nrow=1,ncol=length(sphmc_0.8$omega)))
colnames(new_columns_omega)=name2
result_parameters_0.8=data.frame(task_id=task_id,method='sphmc_0.8',new_columns_beta,new_columns_omega)

new_columns_beta=data.frame(matrix(sphmc_1$beta,nrow=1,ncol=length(sphmc_1$beta)))
colnames(new_columns_beta)=name1
new_columns_omega=data.frame(matrix(sphmc_1$omega,nrow=1,ncol=length(sphmc_1$omega)))
colnames(new_columns_omega)=name2
result_parameters_1=data.frame(task_id=task_id,method='sphmc_1',new_columns_beta,new_columns_omega)

new_columns_beta=data.frame(matrix(sphmc_2$beta,nrow=1,ncol=length(sphmc_2$beta)))
colnames(new_columns_beta)=name1
new_columns_omega=data.frame(matrix(sphmc_2$omega,nrow=1,ncol=length(sphmc_2$omega)))
colnames(new_columns_omega)=name2
result_parameters_2=data.frame(task_id=task_id,method='sphmc_2',new_columns_beta,new_columns_omega)

new_columns_beta=data.frame(matrix(hmc$beta,nrow=1,ncol=length(hmc$beta)))
colnames(new_columns_beta)=name1
new_columns_omega=data.frame(matrix(hmc$omega,nrow=1,ncol=length(hmc$omega)))
colnames(new_columns_omega)=name2
result_parameters_hmc=data.frame(task_id=task_id,method='hmc',new_columns_beta,new_columns_omega)

write.table(result,file='result/water_pca/test_water_pca.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/water_pca/test_water_pca.csv'),append=TRUE)

write.table(result_parameters_0.8,file='result/water_pca/test_water_pca_0.8_parameters.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/water_pca/test_water_pca_0.8_parameters.csv'),append=TRUE)

write.table(result_parameters_1,file='result/water_pca/test_water_pca_1_parameters.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/water_pca/test_water_pca_1_parameters.csv'),append=TRUE)

write.table(result_parameters_2,file='result/water_pca/test_water_pca_2_parameters.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/water_pca/test_water_pca_2_parameters.csv'),append=TRUE)

write.table(result_parameters_hmc,file='result/water_pca/test_water_pca_hmc_parameters.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/water_pca/test_water_pca_hmc_parameters.csv'),append=TRUE)

record_beta=matrix(sphmc_0.8$beta_record,nrow=1,ncol=length(sphmc_0.8$beta_record))
record_0.8_beta=data.frame(task_id=task_id,method='sphmc_0.8',record_beta)
rownames(record_0.8_beta)=name1

file_path=sprintf(
  'result/water_pca/parameter_0.8/beta/test_water_pca_0.8_beta_record_%s.csv',
  task_id
)

write.table(record_0.8_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_0.8$omega_record,nrow=dim(sphmc_0.8$omega_record)[1],ncol=dim(sphmc_0.8$omega_record)[2])
record_0.8_omega=data.frame(task_id=task_id,method='sphmc_0.8',record_omega)
rownames(record_0.8_omega)=name2

file_path=sprintf(
  'result/water_pca/parameter_0.8/omega/test_water_pca_0.8_omega_record_%s.csv',
  task_id
)
write.table(record_0.8_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_0.8$tau_square_record,nrow=1,ncol=length(sphmc_0.8$tau_square_record))
record_0.8_tau_square=data.frame(task_id=task_id,method='sphmc_0.8',record_tau_square)
file_path=sprintf(
  'result/water_pca/parameter_0.8/tau_square/test_water_pca_0.8_tau_square_record_%s.csv',
  task_id
)
write.table(record_0.8_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_0.8$eta_record,nrow=1,ncol=length(sphmc_0.8$eta_record))
record_0.8_eta=data.frame(task_id=task_id,method='sphmc_0.8',record_eta)
file_path=sprintf(
  'result/water_pca/parameter_0.8/eta/test_water_pca_0.8_eta_record_%s.csv',
  task_id
)
write.table(record_0.8_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(sphmc_1$beta_record,nrow=1,ncol=length(sphmc_1$beta_record))
record_1_beta=data.frame(task_id=task_id,method='sphmc_1',record_beta)
rownames(record_1_beta)=name1
file_path=sprintf(
  'result/water_pca/parameter_1/beta/test_water_pca_1_beta_record_%s.csv',
  task_id
)
write.table(record_1_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_1$omega_record,nrow=dim(sphmc_1$omega_record)[1],ncol=dim(sphmc_1$omega_record)[2])
record_1_omega=data.frame(task_id=task_id,method='sphmc_1',record_omega)
rownames(record_1_omega)=name2
file_path=sprintf(
  'result/water_pca/parameter_1/omega/test_water_pca_1_omega_record_%s.csv',
  task_id
)
write.table(record_1_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_1$tau_square_record,nrow=1,ncol=length(sphmc_1$tau_square_record))
record_1_tau_square=data.frame(task_id=task_id,method='sphmc_1',record_tau_square)
file_path=sprintf(
  'result/water_pca/parameter_1/tau_square/test_water_pca_1_tau_square_record_%s.csv',
  task_id
)
write.table(record_1_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_1$eta_record,nrow=1,ncol=length(sphmc_1$eta_record))
record_1_eta=data.frame(task_id=task_id,method='sphmc_1',record_eta)
file_path=sprintf(
  'result/water_pca/parameter_1/eta/test_water_pca_1_eta_record_%s.csv',
  task_id
)
write.table(record_1_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(sphmc_2$beta_record,nrow=1,ncol=length(sphmc_2$beta_record))
record_2_beta=data.frame(task_id=task_id,method='sphmc_2',record_beta)
rownames(record_2_beta)=name1
file_path=sprintf(
  'result/water_pca/parameter_2/beta/test_water_pca_2_beta_record_%s.csv',
  task_id
)
write.table(record_2_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_2$omega_record,nrow=dim(sphmc_2$omega_record)[1],ncol=dim(sphmc_2$omega_record)[2])
record_2_omega=data.frame(task_id=task_id,method='sphmc_2',record_omega)
rownames(record_2_omega)=name2
file_path=sprintf(
  'result/water_pca/parameter_2/omega/test_water_pca_2_omega_record_%s.csv',
  task_id
)
write.table(record_2_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_2$tau_square_record,nrow=1,ncol=length(sphmc_2$tau_square_record))
record_2_tau_square=data.frame(task_id=task_id,method='sphmc_2',record_tau_square)
file_path=sprintf(
  'result/water_pca/parameter_2/tau_square/test_water_pca_2_tau_square_record_%s.csv',
  task_id
)
write.table(record_2_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_2$eta_record,nrow=1,ncol=length(sphmc_2$eta_record))
record_2_eta=data.frame(task_id=task_id,method='sphmc_2',record_eta)
file_path=sprintf(
  'result/water_pca/parameter_2/eta/test_water_pca_2_eta_record_%s.csv',
  task_id
)
write.table(record_2_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(hmc$beta_record,nrow=1,ncol=length(hmc$beta_record))
record_hmc_beta=data.frame(task_id=task_id,method='hmc',record_beta)
rownames(record_hmc_beta)=name1
file_path=sprintf(
  'result/water_pca/parameter_hmc/beta/test_water_pca_hmc_beta_record_%s.csv',
  task_id
)
write.table(record_hmc_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(hmc$omega_record,nrow=dim(hmc$omega_record)[1],ncol=dim(hmc$omega_record)[2])
record_hmc_omega=data.frame(task_id=task_id,method='hmc',record_omega)
rownames(record_hmc_omega)=name2
file_path=sprintf(
  'result/water_pca/parameter_hmc/omega/test_water_pca_hmc_omega_record_%s.csv',
  task_id
)
write.table(record_hmc_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(hmc$tau_square_record,nrow=1,ncol=length(hmc$tau_square_record))
record_hmc_tau_square=data.frame(task_id=task_id,method='hmc',record_tau_square)
file_path=sprintf(
  'result/water_pca/parameter_hmc/tau_square/test_water_pca_hmc_tau_square_record_%s.csv',
  task_id
)
write.table(record_hmc_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(hmc$eta_record,nrow=1,ncol=length(hmc$eta_record))
record_hmc_eta=data.frame(task_id=task_id,method='hmc',record_eta)
file_path=sprintf(
  'result/water_pca/parameter_hmc/eta/test_water_pca_hmc_eta_record_%s.csv',
  task_id
)
write.table(record_hmc_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)




































