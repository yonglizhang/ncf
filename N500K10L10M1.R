


unlink(".RData")
options(scipen=999)

.libPaths("/nas/longleaf/home/yongliz/R/libs/")
library(rugarch)

.libPaths("/nas/longleaf/home/yongliz/R/libs/")
library(fGarch)

.libPaths("/nas/longleaf/home/yongliz/R/libs/")
library(rmgarch)

.libPaths("/nas/longleaf/home/yongliz/R/libs/")
library(earth)

.libPaths("/nas/longleaf/home/yongliz/R/libs/")
library(mgcv)


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

piece.formula <- function(x,y,crossv)  {

                                #k<-1;x<-x[1:n0,];y<-y[1:n0,k]
                                dat<-data.frame(y,x);
                                names(dat)<-c("y",paste("x",1:ncol(x),sep=""));
                                if (crossv==1) {earth.mod<-earth(y~.,data=dat,pmethod="cv",nfold=5,ncross=5)} else {earth.mod<-earth(y~.,data=dat) }
                                
                                
                                zcof1<-strsplit(format(earth.mod, style="pmax"),"\n")[[1]][-1];
                                kterm<-length(zcof1);
                                zcof2<-unlist(strsplit(zcof1,"pmax\\(0\\,"))
                                
                                
                                
                                zcof3<-matrix(unlist(strsplit(unlist(strsplit(trimws(zcof2[2*(1:kterm)]),"\\)"))," - ") ) ,2,kterm);
                                
                                
                                zcof3_sel<-(substr(trimws(zcof3),1,1)=="x");
                                
                                knotss<-as.numeric(as.character(zcof3[!zcof3_sel]));
                                varname<-trimws(zcof3[zcof3_sel]); 
                                
                                dff<-data.frame(varname,knotss);
                                dff<-unique(dff);




                                varsum<-aggregate(knotss~varname,data=dff,length)
                                kkss<-varsum[,2]+1
                                
                                unique.var.name<- varsum$varname
                                
                                K<- length(unique.var.name)

                                FOM<-rep("a",K)
                                for (k in 1:K){

                                              formula.sign <- rep(" - ", varsum[varsum$varname==unique.var.name[k],2] )
                                              formula.sign[dff[dff$varname==unique.var.name[k],2]< 0] <- " + "
                                              FOM[k]<-paste(unique.var.name[k], "+",paste("I(pmax(", unique.var.name[k],formula.sign , abs(dff[dff$varname==unique.var.name[k],2]), ", 0))",collapse = " + ", sep=""))
                                             
                                              }
                                

                                 m<-lm(as.formula(paste("y",paste(FOM,collapse="+"),sep="~")),dat)
                                 
                                 m.res<-resid(m)

                                 slope.coef<-summary(m)$coef[-1,1]
                                 vars<-rep(unique.var.name,kkss)
                                 coef.mt<-data.frame(vars, slope.coef)

                                 
                                 coef.output<-NULL
                                 
                                 for (k in 1:K){

                                               xvar<-unique.var.name[k]
                                               lowb<-c(-Inf,(dff$knotss)[dff$varname==xvar] )
                                               upb <-c((dff$knotss)[dff$varname==unique.var.name[k]] ,Inf)
                                               coef.int<-rep(0,kkss[k])
                                               coef.slp<-(coef.mt$slope.coef)[coef.mt$vars==unique.var.name[k]]
      
                                               
                                               for( jj in 2:kkss[k]){
                                                                    coef.int[jj]<-coef.int[jj-1]+ coef.slp[jj]*(-lowb[jj])
                                                                    coef.slp[jj]<-coef.slp[jj-1]+ coef.slp[jj]
                                                                    }
                                               coef.output<-rbind(coef.output, data.frame(xvar, lowb, upb, coef.int, coef.slp))
                                                
                                               }
                                  output<-list(coef.output, m.res)
                                  output
                                  }
                                          
                                          
                                          
              
gammaf<-function(kkk,xv_seq,xx_seq,yy_seq){


                #xv_seq<-paste("x",1:10,sep=""); xx_seq<-x_cond_mu[n0+n1,]; yy_seq<-x_cond_sigma[n0+n1,]
                L<-length(xv_seq)
                p<-length(kkk)
                gamma_cov<-array(0,c(p,p,L))
                for (f in 1:L){
                
                              #f<-1
                              xv<-xv_seq[f]
                              coef.list<-list(p)
                              lowb_c<-NULL
                              for (ii in 1:p) {
                                              lowb_c<-c(lowb_c, (kkk[[ii]]$lowb)[kkk[[ii]]$xvar==xv][-1])
                                              }
                              lowb_s<-sort(unique(lowb_c))
                              lowb_s<-c(-Inf,lowb_s)
                              upb_s<-c(lowb_s[-1],Inf)
                              
                              
                              
                              KNT<-length(lowb_s)
                              for (ii in 1:p){
                                              y_intercept<-rep(0,KNT)
                                              y_slope<-rep(0,KNT) 
                                              if (sum((kkk[[ii]])$xvar==xv)!=0) {
                                              for (jj in 1:sum((kkk[[ii]])$xvar==xv))
                                                  {
                                                  y_intercept[(lowb_s >=kkk[[ii]]$lowb[kkk[[ii]]$xvar==xv][jj])&(kkk[[ii]]$upb[kkk[[ii]]$xvar==xv][jj]>=upb_s)]<-((kkk[[ii]]$coef.int)[kkk[[ii]]$xvar==xv])[jj]
                                                  y_slope[(lowb_s >=kkk[[ii]]$lowb[kkk[[ii]]$xvar==xv][jj])&(kkk[[ii]]$upb[kkk[[ii]]$xvar==xv][jj]>=upb_s)]<-((kkk[[ii]]$coef.slp)[kkk[[ii]]$xvar==xv])[jj]
                                                  
                                                  } }
                                             coef.list[[ii]]<-data.frame(lowb_s,upb_s,y_intercept,y_slope)
                                             }
                  
                              
                              
                              xx<-xx_seq[f]
                              yy<-yy_seq[f]
                              gamma01_mat<-matrix(0,2*KNT, p)
                              for (ii in 1:p){
                              
                                              gamma0<-rep(0, KNT); gamma1<-rep(0, KNT)
                                              for   (kk in 1:KNT){
                                                                  gamma0[kk]<- (coef.list[[ii]]$y_intercept)[kk]+(coef.list[[ii]]$y_slope)[kk]*xx
                                                                  gamma1[kk]<- (coef.list[[ii]]$y_slope)[kk]*yy
                                                                  }
                                              gamma01<-c(gamma0,gamma1)
                                              gamma01_mat[,ii]<-gamma01
                                              
                                              }
                              #gamma01_mat
                              
                              
                              
                              ###################################################################
                              ###################################################################
                              ###################################################################
                              lbd<-(coef.list[[1]]$lowb_s-xx)/yy; 
                              ubd<-(coef.list[[1]]$upb_s-xx)/yy;
                                        
                              pnorm_ul<-pnorm(ubd)-pnorm(lbd)
                              dnorm_ul<-dnorm(lbd)-dnorm(ubd)
                              
                              lbd_dnorm<-lbd*dnorm(lbd); lbd_dnorm[is.na(lbd_dnorm)]<-0
                              ubd_dnorm<-ubd*dnorm(ubd); ubd_dnorm[is.na(ubd_dnorm)]<-0 
                              d.dnorm_ul<-lbd_dnorm-ubd_dnorm
                              
                                       
                              
                              
                                       
                              var1<-pnorm_ul-pnorm_ul^2
                              var2<-pnorm_ul+ d.dnorm_ul -dnorm_ul^2
                              
                              Theta_cov<-matrix(0,2*KNT, 2*KNT)
                              
                              
                              for (i in 1:(2*KNT)){
                              for (j in 1:(2*KNT)){
                                            if (i!=j &i<= KNT &j<=KNT ) {Theta_cov[i,j]<- -pnorm_ul[i]*pnorm_ul[j]}
                                            if (i!=j &i> KNT &j>KNT )   {Theta_cov[i,j]<- -dnorm_ul[i-KNT]*dnorm_ul[j-KNT]}
                                            
                                            if ((j-i)!=KNT &i<=KNT & j>KNT )  {Theta_cov[i,j]<- -pnorm_ul[i]*dnorm_ul[j-KNT]}
                                            if ((i-j)!=KNT &i> KNT & j<=KNT)  {Theta_cov[i,j]<- -dnorm_ul[i-KNT]*pnorm_ul[j]}
                                            
                                            if ((j-i)==KNT &i<=KNT & j>KNT )  {Theta_cov[i,j]<- dnorm_ul[j-KNT]-pnorm_ul[i]*dnorm_ul[j-KNT]}
                                            if ((i-j)==KNT &i> KNT & j<=KNT)  {Theta_cov[i,j]<- dnorm_ul[i-KNT]-dnorm_ul[i-KNT]*pnorm_ul[j]}
                                            
                                            }
                                            }
                                            
                              diag(Theta_cov)<-c(var1,var2); 
                              #Theta_cov
                              
                              
                              ###################################################################
                              ###################################################################
                              ###################################################################
                              
                              
                  
                              for (ii1 in 1:p){
                              for (ii2 in 1:p){
                  
                                              gamma_cov[ii1,ii2,f]<-t(gamma01_mat[,ii1])%*% Theta_cov %*%(gamma01_mat[,ii2])
                                              
                                              }
                                              }
                                                            
                              
                              }
                  
                              gamma_cov
                              


}

  
################################################################################
################################################################################
################################################################################
################################################################################
###########################Simulation###########################################
################################################################################
################################################################################
corr_true_vec<-NULL
dcc_est_vec<-NULL
ncf_est_vec<-NULL

MMMM<-1
K<-10
L<-10
L0<-L/2
n0<-500
n1<-10

################################################################################
################################################################################

x_cond_sigma<-matrix(0,n0+n1,L)
x_cond_mu<-matrix(0,n0+n1,L)
x_value<-matrix(0,n0+n1,L)


for (f in 1:L){
              
              xspec = garchSpec(model = list(ar=0.5, omega = 1),cond.dist = c("norm")) 
              
              x<-garchSim(xspec, n = n0+n1, extended = TRUE)
                            
              x_cond_sigma[,f]<-x[,"sigma"]
              x_cond_mu[,f]<-x[,"garch"]-x[,"sigma"]*x[,"eps"]
              x_value[,f]<-x[,"garch"]

              }
x<-data.frame(x_value)
names(x)<-paste("X",1:L, sep="")


################################################################################
z_cond_sigma<-matrix(0,n0+n1,L)
z_cond_mu<-matrix(0,n0+n1,L)
z_value<-matrix(0,n0+n1,L)


for (f in 1:L){
              
                zspec = garchSpec(model = list(ar=0.5, omega = 1),cond.dist = c("norm")) 
                
                z<-garchSim(zspec, n = n0+n1, extended = TRUE)
                              
                z_cond_sigma[,f]<-z[,"sigma"]
                z_cond_mu[,f]<-z[,"garch"]-z[,"sigma"]*z[,"eps"]
                z_value[,f]<-z[,"garch"]
  
                }
                
                
z<-data.frame(z_value)
names(z)<-paste("Z",1:L, sep="")




################################################################################
################################################################################


erspec <- garchSpec(model = list(ar=0.5, omega = 0.25),cond.dist = c("norm")) 

er_value<-matrix(0,n0+n1,K)
er_cond_sigma<-matrix(0,n0+n1,K)
er_cond_mu<-matrix(0,n0+n1,K)


y<-matrix(0,n0+n1,K)



a_seq<-matrix(rnorm(K*L)*4,K,L); a_seq[,(L0+1):L]<-0
b_seq<-matrix(rnorm(K*L)*4,K,L); b_seq[,(L0+1):L]<-0
 


for (k in 1:K)  {
                

                er<-garchSim(erspec, n = n0+n1, extended = TRUE)
                
                er_value[,k]<-er[,"garch"]
                er_cond_sigma[,k]<-er[,"sigma"] 
                er_cond_mu[,k]<-er[,"garch"]-er[,"sigma"]*er[,"eps"]
 
                
                ##Linear Model
                if (MMMM==1) {y[,k]<- as.matrix(x)   %*% (a_seq[k,])+er_value[,k]}
                                    
                ##Nonlinear Model
                if (MMMM==2) {y[,k]<- as.matrix(x^2) %*% (a_seq[k,])+er_value[,k]}
                
                
                
                
                }
               
                
Dat<-data.frame(y)
names(Dat)<-paste("Y",1:K,sep="")


Dat0<-Dat[1:n0,]




################################################################################
################################################################################

corr_true<-array(0,c(K,K,n1))

for (jjmm in 1:n1){

covm_true<-array(0, c(K,K,L))
 

for (f in 1:L ){
for (k1 in 1:K){
for (k2 in 1:K){
                #k1<-25;k2<-25;f<-6
                a1<-a_seq[k1,f]
                a2<-a_seq[k2,f]
                b1<-b_seq[k1,f]
                b2<-b_seq[k2,f]
                
                
                if (MMMM==1)  {covm_true[k1,k2,f]<-a1*a2*( x_cond_sigma[n0+jjmm,f]^2)+(k1==k2)*er_cond_sigma[n0+jjmm,k1]*er_cond_sigma[n0+jjmm,k2] } 
 
                if (MMMM==2)  {covm_true[k1,k2,f]<-a1*a2*(4*(x_cond_mu[n0+jjmm,f])^2*(x_cond_sigma[n0+jjmm,f])^2 + 2*x_cond_sigma[n0+jjmm,f]^4)+(k1==k2)*er_cond_sigma[n0+jjmm,k1]*er_cond_sigma[n0+jjmm,k2]  }              
                
                }
                }
                }
cov2cor(apply(covm_true,c(1,2),sum))
corr_true[,,jjmm]
corr_true[,,jjmm]<-cov2cor(apply(covm_true,c(1,2),sum))

}


#b<-rnorm(L); b[1:2]<-0
#y<-rnorm(n0+n1)+as.matrix(x)%*%(b); dat<-data.frame(x,y)
#earth(y~.,data=dat,pmethod="cv",nfold=5); earth(y~.,data=dat ) 



                     
################################################################################
################################################################################
###########################DCC Estimation#######################################
################################################################################
################################################################################




#xspec = ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(garchOrder = c(1,1), model = 'eGARCH'), distribution.model = 'norm')
#uspec1 = multispec(replicate(K, xspec))
#spec1 = dccspec(uspec = uspec1, dccOrder = c(1, 0), distribution = 'mvnorm', robust = T, lag = 1)


#multf <- multifit(uspec1, Dat0)
#dcc.fit.model <- dccfit(spec1, data = Dat0, fit = multf,solver = "solnp")
#dcc.fit.model <- dccfit(spec1, data = Dat0, fit = multf, solver = "solnp")


#xspec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), variance.model = list(garchOrder = c(1,1), model = 'eGARCH'), distribution.model = 'norm')
#uspec.n <- multispec(replicate(K, xspec))
#spec1 <- dccspec(uspec.n, VAR = T, dccOrder = c(1, 1), distribution = 'mvnorm' ,robust = F)
#dcc.fit.model <- dccfit(spec1, data = Dat0)
#dcc_est<-(dccforecast(dcc.fit.model, n.ahead = n1, n.roll = 0, external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster = NULL)@mforecast)$R[[1]]
#dcc_est
dcc_est<-array(NA,c(K,K,n1))

dcc.formula<-function(Dat0) {


uspec<-ugarchspec(mean.model=list(armaOrder = c(1,0)), variance.model=list(garchOrder = c(1,1), model="fGARCH",submodel="GARCH"), distribution.model="norm")
spec1 = dccspec(uspec = multispec( replicate(K, uspec)), VAR=T, lag=1, dccOrder = c(1,1), distribution = "mvnorm",robust = F)
fit1 = dccfit(spec1, data = Dat0)
dcc_est<-(dccforecast(fit1, n.ahead = n1, n.roll = 0, external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster = NULL)@mforecast)$R[[1]]
dcc_est
}
try(dcc_est<-dcc.formula(Dat0),silent=TRUE)





################################################################################
################################################################################
############################NCF Estimation######################################
################################################################################
################################################################################




x0.pred.mean<-matrix(0,L, n1);
x0.pred.sd<-matrix(0,L, n1)
for (f in 1:L){               
              
              x0.fitmodel<-garchFit(formula=~arma(1,0)+garch(1,1), data=x_value[1:n0,f], trace=F, algorithm = c("nlminb"),include.mean=TRUE)
  
              x0.pred<-predict(x0.fitmodel, n.ahead = n1,mse="uncond")
              
              x0.pred.mean[f,]<-x0.pred$meanForecast
              x0.pred.sd[f,]<-x0.pred$standardDeviation
              }


################################################################################
################################################################################
####Crossv==1

ncf1.formula<-function(x,y,x0.pred.mean,x0.pred.sd){
kkk<-list(K); 
y_resid_mat <-matrix(0,n0,K)



for (k in 1:K)   { m_k<-piece.formula(x[1:n0,],y[1:n0,k],1); kkk[[k]]<-m_k[[1]];   y_resid_mat[,k]<-m_k[[2]] }
y1.cond.sigma<-matrix(0,K,n1)
for (k in 1:K)  {
                
                y1fitmodel<-NULL
                y1fitmodel<-garchFit(~arma(1,0)+garch(1,1), data=y_resid_mat[,k],  trace=F, algorithm = c("nlminb"))
                y1.cond.sigma[k,]<-predict(y1fitmodel, n.ahead = n1)$standardDeviation
                
                #y1.garch.resid<-residuals(y1fitmodel,standardize=T)
                #AutocorTest(y1.garch.resid)
                #ArchTest(y1.garch.resid)
                }                
xv_seq<-paste("x",1:L,sep="")

ncf1_est<-array(NA,dim=c(K,K,n1))
for (jjmm in 1:n1){ncf1_est[,,jjmm]<-cov2cor(apply(gammaf(kkk,xv_seq,x0.pred.mean[,jjmm],x0.pred.sd[,jjmm]),c(1,2),sum)+diag(y1.cond.sigma[,jjmm])  )  }
ncf1_est
}

ncf1_est<-array(NA,c(K,K,n1))
try(ncf1_est<-ncf1.formula(x,y,x0.pred.mean,x0.pred.sd),silent=TRUE)


################################################################################
################################################################################
####Crossv==2
ncf2.formula<-function(x,y,x0.pred.mean,x0.pred.sd) {

kkk<-list(K); 
y_resid_mat <-matrix(0,n0,K)



for (k in 1:K)   { m_k<-piece.formula(x[1:n0,],y[1:n0,k],2); kkk[[k]]<-m_k[[1]];   y_resid_mat[,k]<-m_k[[2]] }
y1.cond.sigma<-matrix(0,K,n1)
for (k in 1:K)  {
                
                y1fitmodel<-NULL
                y1fitmodel<-garchFit(~arma(1,0)+garch(1,1), data=y_resid_mat[,k],  trace=F, algorithm = c("nlminb"))
                y1.cond.sigma[k,]<-predict(y1fitmodel, n.ahead = n1)$standardDeviation
                
                #y1.garch.resid<-residuals(y1fitmodel,standardize=T)
                #AutocorTest(y1.garch.resid)
                #ArchTest(y1.garch.resid)
                }                
xv_seq<-paste("x",1:L,sep="")


ncf2_est<-array(NA,dim=c(K,K,n1))
for (jjmm in 1:n1){ ncf2_est[,,jjmm]<-cov2cor(apply(gammaf(kkk,xv_seq,x0.pred.mean[,jjmm],x0.pred.sd[,jjmm]),c(1,2),sum)+diag(y1.cond.sigma[,jjmm])  )  }
ncf2_est
}
ncf2_est<-array(NA,c(K,K,n1))
try(ncf2_est<-ncf2.formula(x,y,x0.pred.mean,x0.pred.sd),silent=TRUE)



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
row_ind<-row(diag(1:K)); col_ind<-col(diag(1:K))
corr_true_vec<-corr_true[row_ind<col_ind]
dcc_est_vec<-(dcc_est)[row_ind<col_ind]
ncf1_est_vec<-(ncf1_est)[row_ind<col_ind]
ncf2_est_vec<-(ncf2_est)[row_ind<col_ind]


output<-cbind(corr_true_vec,dcc_est_vec,ncf1_est_vec,ncf2_est_vec);


print(c(n0, n1, K, L, L0, MMMM, mean(abs(output[,1]-output[,2])), mean(abs(output[,1]-output[,3])),mean(abs(output[,1]-output[,4]))))





writepath<-"/nas/longleaf/home/yongliz/Power/txt3/"
vvv<-round(abs(rnorm(1)*10000000000000))
filename<-paste(vvv, ".n",n0,"K",K, "L",L,"M",MMMM,"txt",sep="")
filename
write.table(output, file=paste(writepath, filename,sep=""), row.names=FALSE,col.names=FALSE)

 
 
 
 
 