f=function(n){
  for(i in c(1:n,(n-1):1)){
    message(noquote(paste(rep(noquote("*"),i),sep='',collapse='')))
  }
}
f(10)
