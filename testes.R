library(MASS)

Box<-function(datos,grupo) 
 {       
 N<-nrow(datos)
 grupo <- as.factor(grupo)
 ng<-length(levels(grupo))
 v<-N-ng 
 p<-length(datos) 
 df <- ((ng - 1) * p * (p + 1))/2 
 
 
 Si<-by(datos,grupo,var) 
 ni<-by(datos,grupo,nrow) 
 ni<-as.vector(ni) 
 ni<-c(unlist(ni)) 
 
 S<-matrix(0,p,p) 
 
 for(j in 1:(length(ni))) 
	 
   {S<-S+ (ni[j]-1)*matrix(unlist(Si[j]),p,p)} 
 
 S<-S/v 
 
 suma<-0 
 
 for(j in 1:(length(ni))) 
	 
   {suma<-suma+ (ni[j]-1)*log(det(matrix(unlist(Si[j]),p,p)))} 
 
 
 Q<- log(det(S))*v-suma 
 
# ajuste: 
 
c1 <- (p * (2 * p + 3) - 1)/(6 * (ng - 1) * (p + 1)) * (sum(1/(ni-1)) - 1/v) 
Q <- c(Box.M = Q, adj.M = (1 -c1) * Q) 
df <- rep(df, 2) 
chi <- data.frame(Statistic = Q,df = df, Pr = 1 -pchisq(Q, df, ncp = 0),row.names = c("Box.M","adj.M")) 
 
print(chi) 
} 



testM<-function(y) 
{n<-nrow(y) 
p<-ncol(y) 
dfchi<-p*(p+1)*(p+2)/6 
q<-diag(n)-(1/n)*matrix(1,n,n) 
y<-as.matrix(y) 
 
s<- t(y) %*% q %*% y 
s<- (1/n)*s 
 
sinv<-ginv(s)    ##sinv<-ginverse(s) ginv existe en MASS 
 
 
gmatriz<-q%*%y%*%sinv%*%t(y)%*%q                   
beta1hat<-(sum(gmatriz*gmatriz*gmatriz))/(n*n) 
beta2hat<-sum(diag(gmatriz*gmatriz))/n 
kappa1<-n*beta1hat/6                                              
kappa2<-(beta2hat-p*(p+2))/sqrt(8*p*(p+2)/n)        
pvalsim<-1-pchisq(kappa1,dfchi) 
pvalkurt<-2*(1-pnorm(abs(kappa2))) 
res<-c("kappa1"=kappa1, "pvalsim"=pvalsim, "kappa2"=kappa2, "pvalkurt"=pvalkurt,"n"=n)
 #"beta1hat"=beta1hat , "beta2hat"=beta2hat,
res 
} 
 
 
 
TestMED<-function(datos,grupo) 
 {     
	grupo <- as.factor(grupo) 
	ng<-length(levels(grupo)) #ng<-length(unlist(factor.names(datos)[col])) 
	p<-length(datos) 
	N<-nrow(datos)        
	v<-N-ng 
 
means<-matrix(0,ng,p)   ##OJO NO hace la mean por variable sino para todo el conj de datos 
for (j in 1:p) 
	{means[,j]<-by(datos[,j],grupo,mean)} 
 
	ni<-by(datos,grupo,nrow)  ##ni<-by(datos[,col],grupos,length) 
	ni<-as.vector(ni)  
	ni<-c(unlist(ni))  
	 
	props <- matrix(ni/sum(ni),ng,1)  
	a<-unlist(means); M<-matrix(a,ng,p)    ###aca taba GRAN parte del problema 
	h <- crossprod(M,props)     
	sca<-scale(M,h,F)  
	SSB <- crossprod(sca*sqrt(ni)) 
 
	X<-datos  
        for(i in 1:N)
            { k <- grupo[i] ;  X[i,]<-X[i,]-M[k,]}
 
	X<-as.matrix(X)  
 
	Si<-by(datos,grupo,var) 
	S<-matrix(0,p,p) 
 
	for(j in 1:(length(ni))) 
		{S<-S+ (ni[j]-1)*matrix(unlist(Si[j]),p,p)} 
			 
	S<-S/v 
	Cov.p <- S 
	di <- sqrt(diag(Cov.p)) 
        Si <- diag(1/di, p)                      

	qx <- svd(X %*% Si) 
 
   		singular.tol = sqrt(.Machine$double.eps)                    
       r <- sum(qx$d > singular.tol*qx$d[1]) 
       scaling.p <- sqrt(v) * (qx$v[,1:r, drop = F] %*% diag(1/qx$d[1:r], r)) 
       if(r < p) scaling.p <- scaling.p %*%t(qx$v[, 1:r, drop =F]) 
       scaling.p  <- Si %*% scaling.p 

        
	esto<-crossprod(t(scaling.p)/sqrt(v)) ; aquello<-SSB 
	lbda <- Re(eigen(esto%*%aquello)$values) 
	lbda <- lbda[abs(lbda) > sqrt(.Machine$double.eps)] 
 
	stat <- c(prod(1/(1 + lbda)), sum(lbda/(1 + lbda)), sum(lbda), lbda[1])

        
	q <- length(lbda) 
	s <- min(p, q) 
	m <- (abs(p - q) - 1)/2 
	n <- (v - p - 1)/2 
	r <- v - (p - q + 1)/2 
	u <- (p * q - 2)/4 
	p2 <- p * p 
	q2 <- q * q 
	pq5 <- p2 + q2 - 5 
	.t <- if(pq5 > 0) sqrt((p2 * q2 - 4)/pq5) else 1 
	lt <- stat[1]^(1/.t) 
	n1 <- 2 * n + s + 1 
	m1 <- 2 * m + s + 1 
	n2 <- 2 * (s * n + 1) 
	sm1 <- s * m1 
	r1 <- max(c(p, q)) 
	r2 <- v - r1 + q 
	v1 <- c(p * q, sm1, sm1, r1) 
	v2 <- c(r * .t - 2 * u, s * n1, n2, r2) 
	.F <- c((v2[1] * (1 - lt))/lt/v1[1], (n1/m1 * stat[2])/(s - stat[2]), (n2 * stat[3])/s/sm1, (stat[4] *r2)/r1) 
 
 
	stats <- data.frame(Statistics = stat, F = .F, df1 = v1, df2 = v2, Pr = 1 - pf(.F, v1, v2),  
	row.names = c("Wilks Lambda", "Pillai Trace", "Hoteling-Lawley Trace", "Roy Greatest Root")) 
 
stats
}                       
    




testes<-function(datos,grupo)
 {grupos<-as.factor(grupo) 
  
 
cat('    ','\n');cat('                  ~~Test de Igualdad de Medias~~ ','\n')
cat('    ','\n')  
  print(TestMED(datos,grupos)) 
cat('    ','\n');cat('    ','\n')


cat('    ','\n');cat('                  ~~Test de Homogeneidad de Variancias~~','\n')
cat('    ','\n')   
 Box(datos, grupos) 
cat('    ','\n');cat('    ','\n')   


cat('    ','\n');cat('                  ~~Test de Mardia de Multinormalidad por grupos~~','\n')  
cat('    ','\n')   
  print(by(datos, grupos,testM)) 
cat('    ','\n');cat('    ','\n')   


} 
   
