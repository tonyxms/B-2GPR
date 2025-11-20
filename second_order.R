#Construct the second-order polynomial basis.
second_order=function(x){
  row=nrow(x)
  col=ncol(x)
  y=matrix(0,nrow=row,ncol=col*(col+1)/2)
  count=1
  for (i in 1:col){
    for (j in i:col){
      y[ ,count]=x[ ,i]*x[ ,j]
      count=count+1
    }
  }
  return (y)
}