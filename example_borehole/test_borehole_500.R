#Example for the Borehole function with original input dimension d = 8, training set size n = 500, and test set size m = 1000.
library(invgamma)
library(lhs)
library(abind)
library(mlegp)
library(laGP)
library(ggplot2)
library(coda)
library(BayesianTools)
library(DiceKriging)
rm(list=ls())
lower<-c(0.05,100,63070,990,63.1,700,1120,9855)
upper<-c(0.15,50000,115600,1110,116,820,1680,12045)

p<-8
m<-1000
n<-500

task_id=as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(task_id+1000)

source('header.R')

## Training data
xx<-maximinLHS(n,p)
minxx<-apply(xx,2,min)
maxxx<-apply(xx,2,max)
xx<-scale(xx,minxx,maxxx-minxx)  ## xx is in the range of [0,1]^p
x<-scale(xx,F,1/(upper-lower))
x<-scale(x,-lower,F)    ## x is in the range of the phyiscal model

y<-2*pi*x[,3]*(x[,4]-x[,6])/(log(x[,2]/x[,1])*(1+2*x[,7]*x[,3]/(log(x[,2]/x[,1])*x[,1]^2*x[,8])+x[,3]/x[,5]))   ## No noise involved. 

dis = 0.06*sd(y)
y=y+runif(n,min=-dis,max=dis)
#y <- y+ rnorm(n, mean = 0, sd=0.02)

## Testing data
x.new<-randomLHS(m,p)
new<-matrix(0,ncol=p,nrow=m)
for (i in 1:m)
  new[i,]<-lower+x.new[i,]*(upper-lower)

truey<-2*pi*new[,3]*(new[,4]-new[,6])/(log(new[,2]/new[,1])*(1+2*new[,7]*new[,3]/(log(new[,2]/new[,1])*new[,1]^2*new[,8])+new[,3]/new[,5]))

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
l=2
rou=0.5

t1=Sys.time()
mlegp_model=mlegp(xx,y,constantMean = 0, min.nugget = 1e-6)
t2=Sys.time()
t_mlegp=difftime(t2,t1,units='mins')
y_pred=predict(mlegp_model,x.new)
rmse_mlegp=sqrt(mean((y_pred-truey)^2))/sd(truey)

eps=sqrt(.Machine$double.eps)
t1=Sys.time()
gpi=newGPsep(xx,y,d=0.1,g=0.1*var(y),dK=TRUE)
mle=mleGPsep(gpi,param="both",tmin=c(eps,eps),tmax=c(10,var(y)))
t2=Sys.time()
t_lagp=difftime(t2,t1,units='mins')
w=predGPsep(gpi,x.new)
rmse_lagp=sqrt(mean((w$mean-truey)^2))/sd(truey)

t1=Sys.time()
kriging_constant=km(~1,design=xx,response=y,covtype="gauss",nugget.estim = TRUE)
t2=Sys.time()
t_kriging_constant=difftime(t2,t1,units='mins')
y_pred=predict(kriging_constant,newdata=x.new,type='UK',checkNames=FALSE)
rmse_kringing_constant=sqrt(mean((y_pred$mean-truey)^2))/sd(truey)

t1=Sys.time()
kriging_linear=km(~.,design=xx,response=y,covtype="gauss",nugget.estim = TRUE)
t2=Sys.time()
t_kriging_linear=difftime(t2,t1,units='mins')
y_pred=predict(kriging_linear,newdata=x.new,type='UK',checkNames=FALSE)
rmse_kringing_linear=sqrt(mean((y_pred$mean-truey)^2))/sd(truey)

t1=Sys.time()
kriging_quadratic=km(~(.)^2+ I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2) + I(X5^2) +I(X6^2)+I(X7^2)+I(X8^2) ,design=xx,response=y,covtype="gauss",nugget.estim = TRUE)
t2=Sys.time()
t_kriging_quadratic=difftime(t2,t1,units='mins')
y_pred=predict(kriging_quadratic,newdata=x.new,type='UK',checkNames=FALSE)
rmse_kringing_quadratic=sqrt(mean((y_pred$mean-truey)^2))/sd(truey)

G=matrix(1,nrow=n,ncol=1)
G=cbind(G,xx)
G=cbind(G,second_order(xx))

newG=matrix(1,nrow=m,ncol=1)
newG=cbind(newG,x.new)
newG=cbind(newG,second_order(x.new))

t1=Sys.time()
sphmc_0.8=bayesian_bridge_sphmc_adp_convergence(xx,y,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=0.8,r)
t2=Sys.time()
t_sphmc_0.8=difftime(t2,t1,units='mins')

sphmc_0.8_pred=bayesian_bridge_predict(xx,y,G,x.new,newG,sphmc_0.8)
rmse_sphmc_0.8 <- sqrt(mean((sphmc_0.8_pred-truey)^2))/sd(truey)

t1=Sys.time()
sphmc_1=bayesian_bridge_sphmc_adp_convergence(xx,y,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=1,r)
t2=Sys.time()
t_sphmc_1=difftime(t2,t1,units='mins')

sphmc_1_pred=bayesian_bridge_predict(xx,y,G,x.new,newG,sphmc_1)
rmse_sphmc_1 <- sqrt(mean((sphmc_1_pred-truey)^2))/sd(truey)

t1=Sys.time()
sphmc_2=bayesian_bridge_sphmc_adp_convergence(xx,y,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=1.8,r)
t2=Sys.time()
t_sphmc_2=difftime(t2,t1,units='mins')

sphmc_2_pred=bayesian_bridge_predict(xx,y,G,x.new,newG,sphmc_2)
rmse_sphmc_2 <- sqrt(mean((sphmc_2_pred-truey)^2))/sd(truey)

t1=Sys.time()
hmc=bayesian_bridge_sphmc_adp_hmc_convergence(xx,y,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,l=l,rou=rou)
t2=Sys.time()
t_hmc=difftime(t2,t1,units='mins')

hmc_pred=bayesian_bridge_predict(xx,y,G,x.new,newG,hmc)
rmse_hmc <- sqrt(mean((hmc_pred-truey)^2))/sd(truey)

result=data.frame(task_id=task_id,rmse_mlegp=rmse_mlegp,t_mlegp=t_mlegp,rmse_lagp=rmse_lagp,t_lagp=t_lagp,rmse_sphmc_0.8=rmse_sphmc_0.8,t_sphmc_0.8=t_sphmc_0.8,rmse_sphmc_1=rmse_sphmc_1,t_sphmc_1=t_sphmc_1,rmse_sphmc_2=rmse_sphmc_2,t_sphmc_2=t_sphmc_2,rmse_hmc=rmse_hmc,t_hmc=t_hmc,rmse_kringing_constant=rmse_kringing_constant,t_kriging_constant=t_kriging_constant,rmse_kringing_linear=rmse_kringing_linear,t_kriging_linear=t_kriging_linear,rmse_kringing_quadratic=rmse_kringing_quadratic,t_kriging_quadratic=t_kriging_quadratic)

name1='Intercept'
for (i in 1:p){
  name1=c(name1,bquote(X[.(i)]))
}
for (i in 1:p){
  for (j in i:p){
    name1=c(name1,bquote(X[.(i)]*X[.(j)]))
  }
}

name2=NULL
for (i in 1:p){
  name2=c(name2,bquote(X[.(i)]))
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

write.table(result,file='result/borehole_500/test_borehole_500.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/borehole_500/test_borehole_500.csv'),append=TRUE)

write.table(result_parameters_0.8,file='result/borehole_500/test_borehole_0.8_parameters_500.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/borehole_500/test_borehole_0.8_parameters_500.csv'),append=TRUE)

write.table(result_parameters_1,file='result/borehole_500/test_borehole_1_parameters_500.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/borehole_500/test_borehole_1_parameters_500.csv'),append=TRUE)

write.table(result_parameters_2,file='result/borehole_500/test_borehole_2_parameters_500.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/borehole_500/test_borehole_2_parameters_500.csv'),append=TRUE)

write.table(result_parameters_hmc,file='result/borehole_500/test_borehole_hmc_parameters_500.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/borehole_500/test_borehole_hmc_parameters_500.csv'),append=TRUE)

record_beta=matrix(sphmc_0.8$beta_record,nrow=dim(sphmc_0.8$beta_record)[1],ncol=dim(sphmc_0.8$beta_record)[2])

record_0.8_beta=data.frame(task_id=task_id,method='sphmc_0.8',record_beta)
rownames(record_0.8_beta)=name1

file_path=sprintf(
  'result/borehole_500/parameter_0.8/beta/test_borehole_0.8_beta_record_500_%s.csv',
  task_id
)

write.table(record_0.8_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_0.8$omega_record,nrow=dim(sphmc_0.8$omega_record)[1],ncol=dim(sphmc_0.8$omega_record)[2])

record_0.8_omega=data.frame(task_id=task_id,method='sphmc_0.8',record_omega)
rownames(record_0.8_omega)=name2

file_path=sprintf(
  'result/borehole_500/parameter_0.8/omega/test_borehole_0.8_omega_record_500_%s.csv',
  task_id
)
write.table(record_0.8_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_0.8$tau_square_record,nrow=1,ncol=length(sphmc_0.8$tau_square_record))
record_0.8_tau_square=data.frame(task_id=task_id,method='sphmc_0.8',record_tau_square)
file_path=sprintf(
  'result/borehole_500/parameter_0.8/tau_square/test_borehole_0.8_tau_square_record_500_%s.csv',
  task_id
)
write.table(record_0.8_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_0.8$eta_record,nrow=1,ncol=length(sphmc_0.8$eta_record))
record_0.8_eta=data.frame(task_id=task_id,method='sphmc_0.8',record_eta)
file_path=sprintf(
  'result/borehole_500/parameter_0.8/eta/test_borehole_0.8_eta_record_500_%s.csv',
  task_id
)
write.table(record_0.8_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(sphmc_1$beta_record,nrow=dim(sphmc_1$beta_record)[1],ncol=dim(sphmc_1$beta_record)[2])
record_1_beta=data.frame(task_id=task_id,method='sphmc_1',record_beta)
rownames(record_1_beta)=name1
file_path=sprintf(
  'result/borehole_500/parameter_1/beta/test_borehole_1_beta_record_500_%s.csv',
  task_id
)
write.table(record_1_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_1$omega_record,nrow=dim(sphmc_1$omega_record)[1],ncol=dim(sphmc_1$omega_record)[2])
record_1_omega=data.frame(task_id=task_id,method='sphmc_1',record_omega)
rownames(record_1_omega)=name2
file_path=sprintf(
  'result/borehole_500/parameter_1/omega/test_borehole_1_omega_record_500_%s.csv',
  task_id
)
write.table(record_1_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_1$tau_square_record,nrow=1,ncol=length(sphmc_1$tau_square_record))
record_1_tau_square=data.frame(task_id=task_id,method='sphmc_1',record_tau_square)
file_path=sprintf(
  'result/borehole_500/parameter_1/tau_square/test_borehole_1_tau_square_record_500_%s.csv',
  task_id
)
write.table(record_1_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_1$eta_record,nrow=1,ncol=length(sphmc_1$eta_record))
record_1_eta=data.frame(task_id=task_id,method='sphmc_1',record_eta)
file_path=sprintf(
  'result/borehole_500/parameter_1/eta/test_borehole_1_eta_record_500_%s.csv',
  task_id
)
write.table(record_1_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(sphmc_2$beta_record,nrow=dim(sphmc_2$beta_record)[1],ncol=dim(sphmc_2$beta_record)[2])
record_2_beta=data.frame(task_id=task_id,method='sphmc_2',record_beta)
rownames(record_2_beta)=name1
file_path=sprintf(
  'result/borehole_500/parameter_2/beta/test_borehole_2_beta_record_500_%s.csv',
  task_id
)
write.table(record_2_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_2$omega_record,nrow=dim(sphmc_2$omega_record)[1],ncol=dim(sphmc_2$omega_record)[2])
record_2_omega=data.frame(task_id=task_id,method='sphmc_2',record_omega)
rownames(record_2_omega)=name2
file_path=sprintf(
  'result/borehole_500/parameter_2/omega/test_borehole_2_omega_record_500_%s.csv',
  task_id
)
write.table(record_2_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_2$tau_square_record,nrow=1,ncol=length(sphmc_2$tau_square_record))
record_2_tau_square=data.frame(task_id=task_id,method='sphmc_2',record_tau_square)
file_path=sprintf(
  'result/borehole_500/parameter_2/tau_square/test_borehole_2_tau_square_record_500_%s.csv',
  task_id
)
write.table(record_2_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_2$eta_record,nrow=1,ncol=length(sphmc_2$eta_record))
record_2_eta=data.frame(task_id=task_id,method='sphmc_2',record_eta)
file_path=sprintf(
  'result/borehole_500/parameter_2/eta/test_borehole_2_eta_record_500_%s.csv',
  task_id
)
write.table(record_2_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(hmc$beta_record,nrow=dim(hmc$beta_record)[1],ncol=dim(hmc$beta_record)[2])
record_hmc_beta=data.frame(task_id=task_id,method='hmc',record_beta)
rownames(record_hmc_beta)=name1
file_path=sprintf(
  'result/borehole_500/parameter_hmc/beta/test_borehole_hmc_beta_record_500_%s.csv',
  task_id
)
write.table(record_hmc_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(hmc$omega_record,nrow=dim(hmc$omega_record)[1],ncol=dim(hmc$omega_record)[2])
record_hmc_omega=data.frame(task_id=task_id,method='hmc',record_omega)
rownames(record_hmc_omega)=name2
file_path=sprintf(
  'result/borehole_500/parameter_hmc/omega/test_borehole_hmc_omega_record_500_%s.csv',
  task_id
)
write.table(record_hmc_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(hmc$tau_square_record,nrow=1,ncol=length(hmc$tau_square_record))
record_hmc_tau_square=data.frame(task_id=task_id,method='hmc',record_tau_square)
file_path=sprintf(
  'result/borehole_500/parameter_hmc/tau_square/test_borehole_hmc_tau_square_record_500_%s.csv',
  task_id
)
write.table(record_hmc_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(hmc$eta_record,nrow=1,ncol=length(hmc$eta_record))
record_hmc_eta=data.frame(task_id=task_id,method='hmc',record_eta)
file_path=sprintf(
  'result/borehole_500/parameter_hmc/eta/test_borehole_hmc_eta_record_500_%s.csv',
  task_id
)
write.table(record_hmc_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)
















