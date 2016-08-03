#options(echo = FALSE)
library(MASS)
library(stats)

source("func.R")

# Import x.data
# Import y.data

args <- commandArgs(trailingOnly = TRUE)

x.data = read.table(args[1])
x.map = x.data[,1]
y.data = read.table(args[2])
y.map = y.data[,1]

#rows are observations, columns are genes(variables)
x.data = t(x.data[,-1])
y.data = t(y.data[,-1])

n.cv <- 5       # Number of cross-validation steps to select sparseness parameters

p <- length(x.map)        # Number of varibles in X
q <- length(y.map)        # Number of varibles in Y

n.sample <- nrow(x.data)  # Sample size (# of observations or chips)

### _______________________________________________________________________________________
### Setting sparseness parameters
### _______________________________________________________________________________________

max.v = 0.4	
max.u = 0.4
step.lambda = 0.01

lambda.v.seq <- seq(0, max.v, by=step.lambda)    # Possible values of sparseness parameters for data Y. Lower bouns should be 0, upper bound can be increased to 0.2.
lambda.u.seq <- seq(0, max.u, by=step.lambda)    # Possible values of sparseness parameters for data X. Lower bouns should be 0, upper bound can be increased to 0.2.

n.lambdas.u <-  length(lambda.u.seq)
n.lambdas.v <-  length(lambda.v.seq)

#lambda.v.matrix <- matrix(rep(lambda.v.seq, n.lambdas.u), nrow=n.lambdas.u, byrow=T)
lambda.v.matrix <- matrix(rep(lambda.v.seq, n.lambdas.u), nrow=n.lambdas.u, byrow=T)
lambda.u.matrix <- matrix(rep(lambda.u.seq, n.lambdas.v), nrow=n.lambdas.u, byrow=F)

ones.p <- rep(1, p)/p
ones.q <- rep(1, q)/q


### _______________________________________________________________________________________
### Analysis
### _______________________________________________________________________________________

#out.file = paste("scca_",n.cv,"_",max.u,"_",max.v,"_",step.lambda,".txt",sep = "")
out.file = "out_scca_CV_gamma.txt"
cat(date(),"\n",file = out.file)
cat("n.cv, max.u, max.v and step.lambda are:\n",n.cv,max.u,max.v,step.lambda,"\n",file = out.file,append = TRUE)
cat("begin select sparseness parameters:\n",file = out.file,append = TRUE)
#cat("begin select sparseness parameters:\n",file = "scca_iteration.txt",append = TRUE)

n.cv.sample <- trunc(n.sample/n.cv)
whole.sample <- seq(1, n.sample)

predict.corr.scca <- matrix(0, nrow=n.lambdas.u, ncol=n.lambdas.v)      # This matrix will contain average test sample correlation for each combination of sparseness parameters

#_______Cross-validation to select optimal combination of sparseness parameters____________
for (i.cv in 1:n.cv)
{
	cat("cross validation",i.cv,":\n",file = out.file,append = TRUE)
        testing.sample <- whole.sample[((i.cv-1)*n.cv.sample+1):(i.cv*n.cv.sample)]
        training.sample <- whole.sample[!whole.sample%in%testing.sample]

        k <- sample.sigma12.function(x.data[training.sample, ], y.data[training.sample, ])

	# Get starting values for singular vectors
        # as column and row means from matrix K
        u.initial <- k %*% ones.q
        u.initial <- u.initial /sqrt(as.numeric(t(u.initial)%*%u.initial))
        v.initial <- t(k) %*% ones.p
        v.initial <- v.initial /sqrt(as.numeric(t(v.initial)%*%v.initial))

        # _______________Data for Predicted correlation (testing sample)_________________

        x.predict <- x.data[testing.sample, ]
        y.predict <- y.data[testing.sample, ]

        # Standardize data
        x.predict <- x.predict - mean(x.predict)
        y.predict <- y.predict - mean(y.predict)

        sigma11.predict <- var(x.predict)
        sigma22.predict <- var(y.predict)

        x.predict <- x.predict %*% diag( 1/sqrt(diag(sigma11.predict)) )
        y.predict <- y.predict %*% diag( 1/sqrt(diag(sigma22.predict)) )

        # ____________Loops for sparseness parameter combinations__________
        for(j.lambda.v in 1:n.lambdas.v)
	{

        	flag.na <- 0

		for(j.lambda.u in 1:n.lambdas.u)
		{
		        lambda.v <- lambda.v.seq[j.lambda.v]    # sparseness parameter for Y
			lambda.u <- lambda.u.seq[j.lambda.u]    # sparseness parameter for X

		        if(flag.na==0)
			{
			        uv <- scca.function(k, u.initial, v.initial, lambda.u, lambda.v)

				vj <- uv$v.new
			        uj <- uv$u.new

			        # Calculate predicted correlation for SCCA
        
			      	predict.corr.scca[j.lambda.u, j.lambda.v] <- predict.corr.scca[j.lambda.u, j.lambda.v] + abs(cor(x.predict%*%uj, y.predict%*%vj))
			        #when either uj or vj or both are zero vector

			        if(is.na(predict.corr.scca[j.lambda.u, j.lambda.v]))
				{
					flag.na <- 1
					cat("NA at",lambda.u.seq[j.lambda.u],"and",lambda.v.seq[j.lambda.v],"\n",file = out.file,append = TRUE)
					#cat(uj,"\n",file = out.file,append = TRUE)
					#cat(vj,"\n",file = out.file,append = TRUE)
				}		
		        }       # close if

			if(flag.na==1)
			{
			        predict.corr.scca[j.lambda.u:n.lambdas.u, j.lambda.v] <- predict.corr.scca[j.lambda.u:n.lambdas.u, j.lambda.v] + NA
			        break
			}
	        }       # close loop on lambda.u
        }       # close loop on lambda.v
}       # close cross-validation loop
# ______________Identify optimal sparseness parameter combination___________

predict.corr.scca[is.na(predict.corr.scca)] <- 0
predict.corr.scca <- predict.corr.scca/n.cv

best.predict.corr.scca <- max(abs(predict.corr.scca), na.rm=T)
best.lambda.v <- lambda.v.matrix[predict.corr.scca==best.predict.corr.scca]
best.lambda.u <- lambda.u.matrix[predict.corr.scca==best.predict.corr.scca]

cat("best average test sample correlation:",best.predict.corr.scca,"is at",best.lambda.u,"and",best.lambda.v,"\n",file = out.file,append = TRUE)

# ______________________________________________________________________________________________________________
# _____Compute singular vectors using the optimal sparseness parameter combination for the whole data___________
# ______________________________________________________________________________________________________________

cat("begin analyze the whole data:","\n",file = out.file,append = TRUE)
#cat("begin analyze the whole data:","\n",file = "scca_iteration.txt",append = TRUE)

k <- sample.sigma12.function(x.data, y.data)

# Get starting values for singular vectors
# as column and row means from matrix K
u.initial <- k %*% ones.q
u.initial <- u.initial /sqrt(as.numeric(t(u.initial)%*%u.initial))
v.initial <- t(k) %*% ones.p
v.initial <- v.initial /sqrt(as.numeric(t(v.initial)%*%v.initial))

uv <- scca.function(k, u.initial, v.initial, best.lambda.u, best.lambda.v)

vj <- uv$v.new  # sparse singular vector (canonical vector for Y)
uj <- uv$u.new  # sparse singular vector (canonical vector for X)

cat("converged after",uv$i,"iterations.","\n",file = out.file,append = TRUE)

corr.scca <- abs(cor(x.data%*%uj, y.data%*%vj)) # canonical correlation for X and Y data

cat("the final overall correlation is",corr.scca,"\n",file = out.file,append = TRUE)
cat("between",sum(uj != 0),"x variables and",sum(vj != 0),"y variables. \n",file = out.file,append = TRUE)
cat("Transcript Factors \n", file = out.file, append = TRUE)

index.u = uj != 0
t1 = uj[index.u]
t2 = x.map[index.u]
index.map = order(t1,decreasing=T)
index.abs.map = order(abs(t1),decreasing=T)
answer.u = data.frame(t1,t2,sort(t1,decreasing=T),t2[index.map],t1[index.abs.map],t2[index.abs.map])
write.table(answer.u,file = out.file,append = TRUE,row.names = TRUE, col.names = FALSE)
cat("\n================================================================================\n",file = out.file,append = TRUE)
cat("\nTarget Genes\n",file = out.file,append = TRUE)

index.v = vj != 0
t1 = vj[index.v]
t2 = y.map[index.v]
index.map = order(t1,decreasing=T)
index.abs.map = order(abs(t1),decreasing=T)
answer.v = data.frame(t1,t2,sort(t1,decreasing=T),t2[index.map],t1[index.abs.map],t2[index.abs.map])
write.table(answer.v,file = out.file,append = TRUE,row.names = TRUE, col.names = FALSE)

