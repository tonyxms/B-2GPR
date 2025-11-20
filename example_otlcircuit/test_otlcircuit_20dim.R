#Example for the OTL Circuit function with increased input dimension d = 20, training set size n = 200, and test set size m = 1000.
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

lower=c(50,25,0.5,1.2,0.25,50)
upper=c(150,70,3,2.5,1.2,300)
p=6
p_fake=20-p
m=1000
n=200

task_id=as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(task_id+1000)

source('header.R')

#train data
xx.train<-maximinLHS(n,p)
minxx.train<-apply(xx.train,2,min)
maxxx.train<-apply(xx.train,2,max)
xx.train<-scale(xx.train,minxx.train,maxxx.train-minxx.train)  
x.train<-scale(xx.train,F,1/(upper-lower))
x.train<-scale(x.train,-lower,F)
y.train=rep(0,n)
for (i in 1:n){
  y.train[i]=otlcircuit(x.train[i, ])
}
dis=0.06*sd(y.train)
y.train=y.train+runif(n,min=-dis,max=dis)
xx.train.fake=maximinLHS(n,p_fake)
xx.train=cbind(xx.train,xx.train.fake)
#test data
x.new<-randomLHS(m,p)
x.test<-matrix(0,ncol=p,nrow=m)
y.test=rep(0,m)
for (i in 1:m){
  x.test[i,]<-lower+x.new[i,]*(upper-lower)
  y.test[i]=otlcircuit(x.test[i, ])
}
x.new.fake=randomLHS(m,p_fake)
x.new=cbind(x.new,x.new.fake)

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
n_adpt=500
r=array(0,2)
r[1]=20
r[2]=5
l=1
rou=0.5

t1=Sys.time()
mlegp_model=mlegp(xx.train,y.train,constantMean = 0, min.nugget = 1e-6)
t2=Sys.time()
t_mlegp=difftime(t2,t1,units='mins')
y_pred=predict(mlegp_model,x.new)
rmse_mlegp=sqrt(mean((y_pred-y.test)^2))/sd(y.test)

eps=sqrt(.Machine$double.eps)
t1=Sys.time()
gpi=newGPsep(xx.train,y.train,d=0.1,g=0.1*var(y.train),dK=TRUE)
mle=mleGPsep(gpi,param="both",tmin=c(eps,eps),tmax=c(10,var(y.train)))
t2=Sys.time()
t_lagp=difftime(t2,t1,units='mins')
w=predGPsep(gpi,x.new)
rmse_lagp=sqrt(mean((w$mean-y.test)^2))/sd(y.test)

t1=Sys.time()
kriging_constant=km(~1,design=xx.train,response=y.train,covtype="gauss",nugget.estim = TRUE)
t2=Sys.time()
t_kriging_constant=difftime(t2,t1,units='mins')
y_pred=predict(kriging_constant,newdata=x.new,type='UK',checkNames=FALSE)
rmse_kringing_constant=sqrt(mean((y_pred$mean-y.test)^2))/sd(y.test)

t1=Sys.time()
kriging_linear=km(~.,design=xx.train,response=y.train,covtype="gauss",nugget.estim = TRUE)
t2=Sys.time()
t_kriging_linear=difftime(t2,t1,units='mins')
y_pred=predict(kriging_linear,newdata=x.new,type='UK',checkNames=FALSE)
rmse_kringing_linear=sqrt(mean((y_pred$mean-y.test)^2))/sd(y.test)

G=matrix(1,nrow=n,ncol=1)
G=cbind(G,xx.train)
#G=cbind(G,second_order(xx.train))

newG=matrix(1,nrow=m,ncol=1)
newG=cbind(newG,x.new)
#newG=cbind(newG,second_order(x.new))

t1=Sys.time()
sphmc_0.8=bayesian_bridge_sphmc_adp_convergence(xx.train,y.train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=0.8,r)
t2=Sys.time()
t_sphmc_0.8=difftime(t2,t1,units='mins')

sphmc_0.8_pred=bayesian_bridge_predict(xx.train,y.train,G,x.new,newG,sphmc_0.8)
rmse_sphmc_0.8 <- sqrt(mean((sphmc_0.8_pred-y.test)^2))/sd(y.test)

t1=Sys.time()
sphmc_1=bayesian_bridge_sphmc_adp_convergence(xx.train,y.train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=1,r)
t2=Sys.time()
t_sphmc_1=difftime(t2,t1,units='mins')

sphmc_1_pred=bayesian_bridge_predict(xx.train,y.train,G,x.new,newG,sphmc_1)
rmse_sphmc_1 <- sqrt(mean((sphmc_1_pred-y.test)^2))/sd(y.test)

t1=Sys.time()
sphmc_2=bayesian_bridge_sphmc_adp_convergence(xx.train,y.train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,q_norm=1.8,r)
t2=Sys.time()
t_sphmc_2=difftime(t2,t1,units='mins')

sphmc_2_pred=bayesian_bridge_predict(xx.train,y.train,G,x.new,newG,sphmc_2)
rmse_sphmc_2 <- sqrt(mean((sphmc_2_pred-y.test)^2))/sd(y.test)

t1=Sys.time()
hmc=bayesian_bridge_sphmc_adp_hmc_convergence(xx.train,y.train,prior_var,G,n_in,n_adpt=n_out*0.8,n_out_min=n_out*1.5,n_out_max=n_out*3,l=l,rou=rou)
t2=Sys.time()
t_hmc=difftime(t2,t1,units='mins')

hmc_pred=bayesian_bridge_predict(xx.train,y.train,G,x.new,newG,hmc)
rmse_hmc <- sqrt(mean((hmc_pred-y.test)^2))/sd(y.test)

result=data.frame(task_id=task_id,rmse_mlegp=rmse_mlegp,t_mlegp=t_mlegp,rmse_lagp=rmse_lagp,t_lagp=t_lagp,rmse_sphmc_0.8=rmse_sphmc_0.8,t_sphmc_0.8=t_sphmc_0.8,rmse_sphmc_1=rmse_sphmc_1,t_sphmc_1=t_sphmc_1,rmse_sphmc_2=rmse_sphmc_2,t_sphmc_2=t_sphmc_2,rmse_hmc=rmse_hmc,t_hmc=t_hmc,rmse_kringing_constant=rmse_kringing_constant,t_kriging_constant=t_kriging_constant,rmse_kringing_linear=rmse_kringing_linear,t_kriging_linear=t_kriging_linear)

name1='Intercept'
for (i in 1:(p+p_fake)){
  name1=c(name1,bquote(X[.(i)]))
}
#for (i in 1:p){
#  for (j in i:p){
#    name1=c(name1,bquote(X[.(i)]*X[.(j)]))
#  }
#}

name2=NULL
for (i in 1:(p+p_fake)){
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


write.table(result,file='result/otlcircuit_200_20/test_otlcircuit_200_20.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/otlcircuit_200_20/test_otlcircuit_200_20.csv'),append=TRUE)

write.table(result_parameters_0.8,file='result/otlcircuit_200_20/test_otlcircuit_0.8_parameters_200_20.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/otlcircuit_200_20/test_otlcircuit_0.8_parameters_200_20.csv'),append=TRUE)

write.table(result_parameters_1,file='result/otlcircuit_200_20/test_otlcircuit_1_parameters_200_20.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/otlcircuit_200_20/test_otlcircuit_1_parameters_200_20.csv'),append=TRUE)

write.table(result_parameters_2,file='result/otlcircuit_200_20/test_otlcircuit_2_parameters_200_20.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/otlcircuit_200_20/test_otlcircuit_2_parameters_200_20.csv'),append=TRUE)

write.table(result_parameters_hmc,file='result/otlcircuit_200_20/test_otlcircuit_hmc_parameters_200_20.csv',sep=',',row.names=FALSE,col.names=!file.exists('result/otlcircuit_200_20/test_otlcircuit_hmc_parameters_200_20.csv'),append=TRUE)

record_beta=matrix(sphmc_0.8$beta_record,nrow=dim(sphmc_0.8$beta_record)[1],ncol=dim(sphmc_0.8$beta_record)[2])

record_0.8_beta=data.frame(task_id=task_id,method='sphmc_0.8',record_beta)
rownames(record_0.8_beta)=name1

file_path=sprintf(
  'result/otlcircuit_200_20/parameter_0.8/beta/test_otlcircuit_0.8_beta_record_200_20_%s.csv',
  task_id
)

write.table(record_0.8_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_0.8$omega_record,nrow=dim(sphmc_0.8$omega_record)[1],ncol=dim(sphmc_0.8$omega_record)[2])
record_0.8_omega=data.frame(task_id=task_id,method='sphmc_0.8',record_omega)
rownames(record_0.8_omega)=name2

file_path=sprintf(
  'result/otlcircuit_200_20/parameter_0.8/omega/test_otlcircuit_0.8_omega_record_200_20_%s.csv',
  task_id
)
write.table(record_0.8_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_0.8$tau_square_record,nrow=1,ncol=length(sphmc_0.8$tau_square_record))
record_0.8_tau_square=data.frame(task_id=task_id,method='sphmc_0.8',record_tau_square)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_0.8/tau_square/test_otlcircuit_0.8_tau_square_record_200_20_%s.csv',
  task_id
)
write.table(record_0.8_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_0.8$eta_record,nrow=1,ncol=length(sphmc_0.8$eta_record))
record_0.8_eta=data.frame(task_id=task_id,method='sphmc_0.8',record_eta)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_0.8/eta/test_otlcircuit_0.8_eta_record_200_20_%s.csv',
  task_id
)
write.table(record_0.8_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(sphmc_1$beta_record,nrow=dim(sphmc_1$beta_record)[1],ncol=dim(sphmc_1$beta_record)[2])
record_1_beta=data.frame(task_id=task_id,method='sphmc_1',record_beta)
rownames(record_1_beta)=name1
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_1/beta/test_otlcircuit_1_beta_record_200_20_%s.csv',
  task_id
)
write.table(record_1_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_1$omega_record,nrow=dim(sphmc_1$omega_record)[1],ncol=dim(sphmc_1$omega_record)[2])
record_1_omega=data.frame(task_id=task_id,method='sphmc_1',record_omega)
rownames(record_1_omega)=name2
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_1/omega/test_otlcircuit_1_omega_record_200_20_%s.csv',
  task_id
)
write.table(record_1_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_1$tau_square_record,nrow=1,ncol=length(sphmc_1$tau_square_record))
record_1_tau_square=data.frame(task_id=task_id,method='sphmc_1',record_tau_square)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_1/tau_square/test_otlcircuit_1_tau_square_record_200_20_%s.csv',
  task_id
)
write.table(record_1_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_1$eta_record,nrow=1,ncol=length(sphmc_1$eta_record))
record_1_eta=data.frame(task_id=task_id,method='sphmc_1',record_eta)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_1/eta/test_otlcircuit_1_eta_record_200_20_%s.csv',
  task_id
)
write.table(record_1_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(sphmc_2$beta_record,nrow=dim(sphmc_2$beta_record)[1],ncol=dim(sphmc_2$beta_record)[2])
record_2_beta=data.frame(task_id=task_id,method='sphmc_2',record_beta)
rownames(record_2_beta)=name1
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_2/beta/test_otlcircuit_2_beta_record_200_20_%s.csv',
  task_id
)
write.table(record_2_beta,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(sphmc_2$omega_record,nrow=dim(sphmc_2$omega_record)[1],ncol=dim(sphmc_2$omega_record)[2])
record_2_omega=data.frame(task_id=task_id,method='sphmc_2',record_omega)
rownames(record_2_omega)=name2
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_2/omega/test_otlcircuit_2_omega_record_200_20_%s.csv',
  task_id
)
write.table(record_2_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(sphmc_2$tau_square_record,nrow=1,ncol=length(sphmc_2$tau_square_record))
record_2_tau_square=data.frame(task_id=task_id,method='sphmc_2',record_tau_square)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_2/tau_square/test_otlcircuit_2_tau_square_record_200_20_%s.csv',
  task_id
)
write.table(record_2_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(sphmc_2$eta_record,nrow=1,ncol=length(sphmc_2$eta_record))
record_2_eta=data.frame(task_id=task_id,method='sphmc_2',record_eta)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_2/eta/test_otlcircuit_2_eta_record_200_20_%s.csv',
  task_id
)
write.table(record_2_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_beta=matrix(hmc$beta_record,nrow=dim(hmc$beta_record)[1],ncol=dim(hmc$beta_record)[2])
record_hmc_beta=data.frame(task_id=task_id,method='hmc',record_beta)
rownames(record_hmc_beta)=name1
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_hmc/beta/test_otlcircuit_hmc_beta_record_200_20_%s.csv',
  task_id
)
write.table(record_hmc_beta,file=file_path,row.names=TRUE,col.names=FALSE,append=TRUE)

record_omega=matrix(hmc$omega_record,nrow=dim(hmc$omega_record)[1],ncol=dim(hmc$omega_record)[2])
record_hmc_omega=data.frame(task_id=task_id,method='hmc',record_omega)
rownames(record_hmc_omega)=name2
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_hmc/omega/test_otlcircuit_hmc_omega_record_200_20_%s.csv',
  task_id
)
write.table(record_hmc_omega,file=file_path,sep=',',row.names=TRUE,col.names=FALSE,append=TRUE)

record_tau_square=matrix(hmc$tau_square_record,nrow=1,ncol=length(hmc$tau_square_record))
record_hmc_tau_square=data.frame(task_id=task_id,method='hmc',record_tau_square)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_hmc/tau_square/test_otlcircuit_hmc_tau_square_record_200_20_%s.csv',
  task_id
)
write.table(record_hmc_tau_square,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)

record_eta=matrix(hmc$eta_record,nrow=1,ncol=length(hmc$eta_record))
record_hmc_eta=data.frame(task_id=task_id,method='hmc',record_eta)
file_path=sprintf(
  'result/otlcircuit_200_20/parameter_hmc/eta/test_otlcircuit_hmc_eta_record_200_20_%s.csv',
  task_id
)
write.table(record_hmc_eta,file=file_path,sep=',',row.names=FALSE,col.names=FALSE,append=TRUE)




















