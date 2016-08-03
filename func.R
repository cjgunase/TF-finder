# This function computes sparse version of the first singular vectors of matrix K

# The number of selected variables, i.e. variables with non-zero entries
# in computed singular vectors, is controlled by the sparseness parameters.
# Increasing the sparseness parameter will decrease the number of selected
# variables.


# Parameters:
# k is p by q covariance matrix for standardized data sets X and Y
# u.initial p by 1 vector of starting values for left singular vector
# v.initial q by 1 vector of starting values for right singular vector
# lambda.u is the sparseness parameter for left singular vector
# lambda.v is the sparseness parameter for right singular vector

scca.function <- function(k, u.initial, v.initial, lambda.u, lambda.v)
{
        i <- 0          # number of iterations used by SCCA
        eps <- 0.001    # convergence criterion
        max.iter <- 50  # maximum number of iterations
        diff.u <- eps*10
        diff.v <- eps*10

	#cat("begin iteration:\n",file = "scca_iteration.txt",append=TRUE)

        while ((i < max.iter) & ((diff.u > eps) || (diff.v > eps)) )
        {
        	i <- i+1

		#cat("iteration",i,":\n",file = "scca_iteration.txt",append=TRUE)

	        # Update left singular vector

        	vx <-  k %*% v.initial
	        length.vx <- as.numeric(sqrt(t(vx)%*%vx))
        	if(length.vx==0) length.vx <- 1
	        vx <- vx / length.vx
        	u.new <- abs(vx) - 0.5*lambda.u
	        u.new <- (u.new + abs(u.new))/2
        	u.new <- u.new*sign(vx)
	        length.u.new <- as.numeric(sqrt(t(u.new)%*%u.new))
        	if(length.u.new==0) length.u.new <- 1
	        u.new <- u.new / length.u.new

	
        	# Update right singular vector

	        ux <- t(k) %*% u.new
        	length.ux <- as.numeric(sqrt(t(ux)%*%ux))
	        if(length.ux==0) length.ux <- 1
        	ux <- ux / length.ux
	        v.new <- abs(ux) - 0.5*lambda.v
        	v.new <- (v.new + abs(v.new))/2
	        v.new <- v.new * sign(ux)
        	length.v.new <- as.numeric(sqrt(t(v.new)%*%v.new))
	        if(length.v.new==0) length.v.new <- 1
        	v.new <- v.new / length.v.new
	        # Convergence measures

        	diff.v <- max(abs(v.initial - v.new))
	        diff.u <- max(abs(u.initial - u.new))

		#cat("u is",u.new,"\n",file = "scca_iteration.txt",append=TRUE)
		#cat("v is",v.new,"\n",file = "scca_iteration.txt",append=TRUE)
		#cat("Difference is",max(diff.u,diff.v),"\n",file = "scca_iteration.txt",append=TRUE)

        	v.initial <- v.new
	        u.initial <- u.new
        }

	# Report the results:
	# u.new is computed left singular vector
	# v.new is computed right singular vector
	# i is the number of iterations used by SCCA

	return(list(u.new=u.new, v.new=v.new, i=i))
}

adaptive.scca.function <- function(k, u.initial, v.initial, lambda.u, lambda.v, u.svd, v.svd, gamma)
{
        i <- 0          # number of iterations used by SCCA
        eps <- 0.001    # convergence criterion
        max.iter <- 50  # maximum nuber of iterations
        diff.u <- eps*10
        diff.v <- eps*10

        # Adjusting weights in for soft-thresholding
        lambda.u.svd <- 1/(abs(u.svd)^gamma)
        lambda.u.svd <- lambda.u.svd*lambda.u

        lambda.v.svd <- 1/(abs(v.svd)^gamma)
        lambda.v.svd <- lambda.v.svd*lambda.v

        while ((i < max.iter) & ((diff.u > eps) || (diff.v > eps)) )
        {
                i <- i+1

                # Update left singular vector

                vx <-  k %*% v.initial
                length.vx <- as.numeric(sqrt(t(vx)%*%vx))
                if(length.vx==0) length.vx <- 1
                vx <- vx / length.vx
                u.new <- abs(vx) - 0.5*lambda.u.svd
                u.new <- (u.new + abs(u.new))/2
                u.new <- u.new*sign(vx)
                length.u.new <- as.numeric(sqrt(t(u.new)%*%u.new))
                if(length.u.new==0) length.u.new <- 1
                u.new <- u.new / length.u.new

                # Update right singular vector

                ux <- t(k) %*% u.new
                length.ux <- as.numeric(sqrt(t(ux)%*%ux))
                if(length.ux==0) length.ux <- 1
                ux <- ux / length.ux
                v.new <- abs(ux) - 0.5*lambda.v.svd
                v.new <- (v.new + abs(v.new))/2
                v.new <- v.new * sign(ux)
                length.v.new <- as.numeric(sqrt(t(v.new)%*%v.new))
                if(length.v.new==0) length.v.new <- 1
                v.new <- v.new / length.v.new

                # Convergence measures

                diff.v <- max(abs(v.initial - v.new))
                diff.u <- max(abs(u.initial - u.new))

                v.initial <- v.new
                u.initial <- u.new
        }

# Report the results:
# u.new is computed left singular vector
# v.new is computed right singular vector
# i is the number of iterations used by adaptive SCCA

	return(list(u.new=u.new, v.new=v.new, i=i))
}

############  SOURCE FOR COVARIANCE MATRIX FUNCTION   #############################
# see source("sample_cov_function.R") above
###################################################################################
# Calculating sample covariance function

sample.sigma12.function <- function(x, y)
{

	#centerize x and y
	#wondering why they are the same result,though different ways?
	#n.x = nrow(x) #x and y should have the same # of rows
	#x <- x-matrix(rep(colMeans(x),n.x),nrow=n.x,byrow=T)
	#y <- y-matrix(rep(colMeans(y),n.x),nrow=n.x,byrow=T)

	# standardize data
	x <- x - mean(x)
	y <- y - mean(y)

	# Sample variance-covariance matrices
	sigma11 <- var(x)
	sigma22 <- var(y)

	x <- x %*% diag( 1/sqrt(diag(sigma11)) )
	y <- y %*% diag( 1/sqrt(diag(sigma22)) )

	return(cov(x,y))
}

