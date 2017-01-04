#### get data mu for each gene #####
get_muhat <- function(y, offset){
  	return(apply(y/exp(offset), 1, mean))
}

### get sigma ###
getsigma <- function(mu, mu_glm){
  	return(sd(log(mu)-log(mu_glm)))
}


#### MLE to estimate genewise dispersion ###
get_genewise_dispersion_mle_mu <- function(y, mu, offset, cores){
	ngene <- nrow(y)
	nsam <- ncol(y)
	muA <- exp(offset)*mu
	tmp <- cbind(y, muA)
  	getdispersionMLE <- function(data){
    #print(1)
    		n <- length(data)
    		y <- data[1:(n/2)]
    		mu <- data[(n/2+1):n]
   		obj <- function(phi){
      			-sum( lgamma(y + 1/phi) - lgamma(1/phi) -lgamma(y+1) - 1/phi*log(1+mu*phi) + y*( log(mu) - log(1/phi+mu) ) )  ###optimize the probability mass function
    		}
    		return(optimize(obj, interval=c(10^(-20), 10^20))$minimum)
  	}
	x <- lapply(apply(tmp, 1, FUN=list), unlist) 
	dispersion <- unlist(mclapply(x, getdispersionMLE ,mc.cores=cores))
  	#dispersion <- apply(tmp, 1, getdispersionMLE)
  	return(dispersion)
}


### other
get_binned_dispersion_mle_mu <- function(y, mu, b, offset, nbin){
  ngene <- dim(y)[1]
  nsam <- dim(y)[2]
  n <- (ngene*nsam-length(b))/nbin
  muA <- exp(offset)*mu

  order_mu <- order(mu, decreasing=T)
  bin_num <- floor(ngene/nbin)
  labs_order <- c( rep(1:nbin,each=bin_num), rep(nbin,ngene-bin_num*nbin) )
  labs <- rep(NA, ngene)
  
  ###rearrange y, muA
  y_R <- c()
  muA_R <- c()
  for ( bin in 1:(nbin-1) ){
    y_R <- rbind(y_R, as.numeric(y[order_mu[which(labs_order==bin)],]) )
    muA_R <- rbind(muA_R, as.numeric(muA[order_mu[which(labs_order==bin)],]) )
    labs[order_mu[which(labs_order==bin)]] <- bin
  }
  bin <- nbin
  y_R_last <- as.numeric(y[order_mu[which(labs_order==bin)],])
  muA_R_last <- as.numeric(muA[order_mu[which(labs_order==bin)],])
  labs[order_mu[which(labs_order==bin)]] <- bin
  
  getdispersionMLE <- function(data){
    #print(1)
    n <- length(data)
    y <- data[1:(n/2)]
    mu <- data[(n/2+1):n]
    obj <- function(phi){
      -sum( lgamma(y+1/phi) - lgamma(1/phi) - 1/phi*log(1+mu*phi) + y*( log(mu) - log(1/phi+mu) ) )
    }
    return(optimize(obj, interval=c(10^(-20), 10^20))$minimum)
  }
  tmp <- cbind(y_R, muA_R)
  dispersion <- apply(tmp, 1, getdispersionMLE)
  dispersion_last <- getdispersionMLE(c(y_R_last, muA_R_last))
  
  dispersion <- c(dispersion, dispersion_last)
  dispersionA <- rep(NA, ngene)
  for (bin in 1:nbin){
    dispersionA[which(labs==bin)] <- dispersion[bin]
  }
  return(dispersionA)  
}



### get optimized mu ########
# mu1 = mu_hat
# mu2 = mu_hat_glm
getmu_post_optimize <- function(y, offset, mu1, mu2, phi, sigma, span,cores){
	ngene <- nrow(y)
	nsam <- ncol(y)  
  	if (length(phi)==1){
    		phi <- rep(phi, ngene)
  	}
  	tmp <- cbind(y, offset, mu1, mu2, phi)
  	getmu <- function(data, sigma, span){
    	#print(1)
    		n <- length(data)
    		disp <- data[n]
    		mu2 <- data[n-1]
    		mu1 <- data[n-2]
    		y <- data[1:((n-3)/2)]
    		offset <- exp(data[((n-3)/2+1):(n-3)])
    		obj <- function(mu){
      			if (disp==0){
        	       		return((log(mu)-log(mu2))^2/(2 * sigma^2) - sum( y*log(mu * offset)-offset*mu -lgamma(y+1)) )    
      			} else {
      				return ((log(mu)-log(mu2))^2/(2 * sigma^2) -sum( lgamma(y + 1/disp) - lgamma(1/disp) -lgamma(y+1) - 1/disp*log(1+mu*disp*offset) + y*( log(mu*offset) - log(1/disp+mu*offset) ) ) )	
			}
    		}
    		return(optimize(obj, interval=c(10^(-10), max(mu1, mu2)*span))$minimum)
  	}
	x <- lapply(apply(tmp, 1, FUN=list), unlist)
        mu <- unlist(mclapply(x, function(x) getmu(x, sigma,span) ,mc.cores=cores))
  	#mu <- apply(tmp, 1, getmu, sigma, span)
  	return(mu)
}

###### from phi estimate mu
get_mu_hat_mle <- function(y, phi, offset,cores){
	ngene <- nrow(y)
        nsam <- ncol(y)
        tmp <- cbind(y,offset,phi)

        getmuMLE <- function(data){
                n <- length(data)-1
                y <- data[1:(n/2)]
                offset <- data[(n/2+1):n]
		phi.g = data[(n+1)]

                obj <- function(x){
			 mu <- exp(offset)*x
                        -sum( lgamma(y + 1/phi.g) - lgamma(1/phi.g) -lgamma(y+1) - 1/phi.g*log(1+mu*phi.g) + y*( log(mu) - log(1/phi.g+mu) ) )  ###optimize the probability mass function
                }
                return(optimize(obj, interval=c(10^(-20), 10^20))$minimum)
        }
        
	x <- lapply(apply(tmp, 1, FUN=list), unlist)
        mu <- unlist(mclapply(x, getmuMLE ,mc.cores=cores))

	return(mu)
	
}

####### apply NB for mu & phi #####
get_mu_phi <- function(y, offset){
	nsam <- ncol(y)
	tmp = cbind(y,offset)
	get_nb <- function(tmp){
		data = data.frame(y=tmp[1:nsam], exp=rep(0,nsam), offset=tmp[(nsam+1):length(tmp)])
		fit = glm.nb(y~exp+offset(offset),data=data)
		return(c(exp(fit$coefficients[1]), 1/fit$theta))
	}
	mu_phi <- apply(tmp,1,get_nb)
	return(mu_phi)
}




###############
### get optimized mu ########
# mu1 = mu_hat
# mu2 = mu_hat_glm
getmu_post_optimize_sample <- function(y, offset, mu1, mu2, phi, sigma, span, cores){
        ngene <- nrow(y)
        nsam <- ncol(y)
        if (length(phi)==1){
                phi <- rep(phi, ngene)
        }
        tmp <- cbind(y, offset, mu1, mu2, phi)
        getmu <- function(data, sigma, span){
        #print(1)
                n <- length(data)
                disp <- data[n]
                mu2 <- data[n-1]
                mu1 <- data[n-2]
                y <- data[1:((n-3)/2)]
                offset <- data[((n-3)/2+1):(n-3)]
		non.zeros = which(offset != 0)
		if (length(non.zeros) ==0){
			return(mu2)
		}else{
                	obj <- function(mu){
                        	if (disp==0){
                                	return((log(mu)-log(mu2))^2/(2 * sigma^2) - sum( y[non.zeros]*log(mu * offset[non.zeros])-offset[non.zeros]*mu -lgamma(y[non.zeros]+1)) )
                        	} else {
                                	return ((log(mu)-log(mu2))^2/(2 * sigma^2) -sum( lgamma(y[non.zeros] + 1/disp) - lgamma(1/disp) -lgamma(y[non.zeros]+1) - 1/disp*log(1+mu*disp*offset[non.zeros]) + y[non.zeros]*( log(mu*offset[non.zeros]) - log(1/disp+mu*offset[non.zeros]) ) ) )
                        	}
                	}
                	return(optimize(obj, interval=c(10^(-10), max(mu1, mu2)*span))$minimum)
		}
        }
	
	x <- lapply(apply(tmp, 1, FUN=list), unlist)
        mu <- unlist(mclapply(x, function(x) getmu(x, sigma,span) ,mc.cores=cores))

        #mu <- apply(tmp, 1, getmu, sigma, span)
        return(mu)
}

#### get data mu for each gene #####
get_muhat_sample <- function(y, offset){
        return(apply(y/offset, 1, function(x) mean(x, na.rm=T)))
}



#### MLE to estimate genewise dispersion ###
get_genewise_dispersion_mle_mu_sample<- function(y, mu, offset, cores){
        tmp <- cbind(y, offset, mu)
        
	getdispersionMLE <- function(data){
    #print(1)
                n <- length(data)
                y <- data[1:((n-1)/2)]
                offset <- data[((n-1)/2+1):(n-1)]
		m = data[n]

		non.zeros = which(offset !=0)
		y = y[non.zeros]
		mu = offset[non.zeros] * m

                obj <- function(phi){
                        -sum( lgamma(y + 1/phi) - lgamma(1/phi) -lgamma(y+1) - 1/phi*log(1+mu*phi) + y*( log(mu) - log(1/phi+mu) ) )  ###optimize the probability mass function
                }
                return(optimize(obj, interval=c(10^(-20), 10^20))$minimum)
        }
        x <- lapply(apply(tmp, 1, FUN=list), unlist)
        dispersion <- unlist(mclapply(x, getdispersionMLE ,mc.cores=cores))
        #dispersion <- apply(tmp, 1, getdispersionMLE)
        return(dispersion)
}


###### from phi estimate mu
get_mu_hat_mle_sample <- function(y, phi, offset,cores){
        tmp <- cbind(y,offset,phi)

        getmuMLE <- function(data){
                n <- length(data)-1
                y <- data[1:(n/2)]
                offset <- data[(n/2+1):n]
        	non.zeros = which(offset !=0)
		y = y[non.zeros]
		offset = offset[non.zeros]

	        phi.g = data[(n+1)]

                obj <- function(x){
                         mu <- offset *x
                        -sum( lgamma(y + 1/phi.g) - lgamma(1/phi.g) -lgamma(y+1) - 1/phi.g*log(1+mu*phi.g) + y*( log(mu) - log(1/phi.g+mu) ) )  ###optimize the probability mass function
                }
                return(optimize(obj, interval=c(10^(-20), 10^20))$minimum)
        }
	x <- lapply(apply(tmp, 1, FUN=list), unlist)
	mu <- unlist(mclapply(x, getmuMLE ,mc.cores=cores))
        #mu <- apply(tmp, 1, getmuMLE)
        return(mu)

}
