suppressMessages(library(abind))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(parallel))
suppressMessages(library("data.table"))
suppressMessages(library(poibin))
suppressMessages(library(survcomp))
suppressMessages(library(methods))
suppressMessages(library(MASS))
source("Code/functions.R")


## check arguments ##
args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  	stop("Usage: Rscript **.R input_maf_file output_file f\n input_maf_file: file containing all somatic SNVs and Indels\n output_file: \n f: fraction of genes to be considered in constructing background model\n", call.=FALSE)
} else if (length(args) > 0) {
	maf.file <- args[1]
	fraction <- as.numeric(args[2])
	output <- args[3]
	method <- args[4]
	cores <- as.numeric(args[5])
	beds <- args[6]	#### these are bases for samples	
	og.num = as.numeric(args[7])
	bin.size =as.numeric(args[8])
	interval.size =as.numeric(args[9])
	pathway.file =args[10]
}






## required files ##
exome.file="/net/kodiak/volumes/delta/shared/home/yao_new/MADGiC_mod/data/exome_hg19_vep_Feb11_2016.Rdata"   #### hg19 exome file
covar.file="/net/kodiak/volumes/delta/shared/home/yao_new/MADGiC_mod/data/gene.covar.txt"   ### gene properties
mutation.context = "/net/kodiak/volumes/delta/shared/home/yao_new/MADGiC_mod/data/mutation_context_96.txt"   ### 96 mutation contexts
exome.gene="/net/kodiak/volumes/delta/shared/home/yao_new/MADGiC_mod/data/exome_gene_vep_Feb11.Rdata"   ### exome info based on gene levels 

### cancer gene census (only for annotation purpose) ###
cgc <- read.table("/net/kodiak/volumes/river/shared/users/yao/Database/cancer_gene_census.txt",header=T,sep="\t")



##############################
########## functions #########
##############################

# get mutations on exomes, get ride of bad annotations #
get.mutation.exome <- function(exome,mutab){
        res = NULL
        mutab[,2] = gsub("X",23,mutab[,2])
        mutab[,2] = gsub("Y",24,mutab[,2])
        mutab = mutab[mutab[,2] %in% 1:24,]
        for (i in 1:23){
                X = data.table(exome[[i]])
                setkey(X,V1,V7)
                sel.mutab = mutab[mutab[,2]==i,]
                which = !is.na(X[J(sel.mutab$Start_Position,sel.mutab$Ensembl_gene_id),]$V2)
                res = rbind(res,sel.mutab[which,])
        }
        return(res)
}

# get background count for silent / non-silent mutations #
preprocess.BM<-function(X,gene, mut.context)
{
	type.const = rbind(rep(0,96,),rep(0,96))
	colnames(type.const) = 1:96
	sel = data.table(X[X[,1] %in% gene,])
	sel$ind = as.numeric(sel$ind >0)
	res = sel[, sum(V1), by = 'tag,ind']
	
	sil = res[res$ind ==0,]
	nonsil = res[res$ind == 1,]

	e = sil[match(colnames(type.const),sil$tag),]$V1
	f = nonsil[match(colnames(type.const),nonsil$tag),]$V1
 	e[is.na(e)] = 0
	f[is.na(f)] = 0
	
	type.const = type.const + rbind(e,f)
  	return(type.const)
}

# get count for silent/non-silent mutations in samples#
mut.type.converter<- function(gene,mutab,mut.context)
{
        res = rep(0,96)
        mutab.snp= mutab[ mutab[,"Context"] != "INDEL" & mutab[,"Ensembl_gene_id"] %in% gene ,]

        mutab.snp[,"Context"] = mut.context[match(mutab.snp[,"Context"],mut.context[,3]),5]
        x = table(mutab.snp[,"Context"])
        y = x[match(1:96,names(x))]
        y[is.na(y)] = 0
        names(y) = 1:96
        res = res + y
  	return(res)
}

# calculate the relative mutation frequency for each tri-nucleotide context #
fit.background.model<-function(mutab, nonsil.mut.type.sampl.sum, sil.mut.type.sampl.sum, nonsil.type.const, sil.type.const,mut.context){
	x = which(sil.type.const[1,] >= 1000 & sil.mut.type.sampl.sum > 0)
	r <-  mean( (nonsil.mut.type.sampl.sum[x]/ nonsil.type.const[2,x]) * (sil.type.const[1,x] /sil.mut.type.sampl.sum[x]))      ######### relative to silent mutations 
	p.all <- ( nonsil.mut.type.sampl.sum + sil.mut.type.sampl.sum ) / ( sil.type.const[1,] + r * nonsil.type.const[2,] )  # part of phat

	ref <- p.all[1]  # reference cate is type 1
	p <- rep(0,96)
     	for (m in 1:96) {
        	p[m] <- p.all[m]/ref
     	}

  ## calculate p_{indel}, which will be used later for calculation of p_{frameshift} and p_{inframe}
	if(sum(  mutab[,4]=="In_frame" | mutab[,4]=="Frame_shift"  )>0){
        	mutab.indel=   mutab[ mutab[,4]=="In_frame" | mutab[,4]=="Frame_shift",]
        	pos = (sum(nonsil.type.const[1,]) + sum(nonsil.type.const[2,]))/3
        	ind <-  nrow(mutab.indel)/pos/r/ref  #### p_{indel}
    	}else{ind <- 0}
    	p=c(p,ind)
    	names(p) = c(1:96,"indel")
    	return(list(p,r,ref))
}


# fit background model with negative-binomial regression #
fit.model <- function(exome.gene,nonsil.mutab,sil.mutab,sample.name,nosam,gene.rep.expr,gid,p,S, which, method, cores){
        X = data.table(data.frame(gene = exome.gene[,1], p=p[as.numeric(exome.gene[,2])]* as.numeric(exome.gene[,4]), ind = as.numeric(exome.gene[,3]>0)))
        gene.prob = X[,sum(p), by = 'gene,ind']
        setkey(gene.prob,gene,ind)
	
        # sample count # 
	mut.count = aggregate(Start_Position ~ Tumor_Sample_Barcode, rbind(sil.mutab,nonsil.mutab),length)
        rownames(mut.count) = mut.count[,1]


        # for all genes with silent mutations 
        sel.gene = intersect(rownames(gene.rep.expr), sil.mutab[,1])
	offset = log(as.matrix(as.numeric(gene.prob[J(sel.gene,0),]$V1)) %*% t(as.matrix(as.numeric(mut.count[sample.name,2]))))   #### offset in log scale 
        count = matrix(0,nrow = length(gid), ncol = S)
        rownames(count) = gid
        colnames(count) = sample.name
        sil.count = aggregate(Start_Position ~ Ensembl_gene_id + Tumor_Sample_Barcode, sil.mutab,length)
        for (i in 1:nrow(sil.count)){
                count[as.character(sil.count[i,1]),as.character(sil.count[i,2])] = sil.count[i,3]
        }
        y = unname(count[sel.gene,])   #### y 
        m= data.frame(epc = as.numeric(gene.prob[J(sel.gene,0),]$V1), exp = as.numeric(gene.rep.expr[sel.gene,3]), rep = as.numeric(gene.rep.expr[sel.gene,4]), hic = as.numeric(gene.rep.expr[sel.gene,5]), or = factor(gene.rep.expr[sel.gene,6], levels =1:2))
        design = model.matrix(~ exp + rep + hic + or, m) ### design
	
	o_scale <- max(offset)
	offsetS <- offset - o_scale ### rescale 
	
	ysum <- apply(y, 1, sum)
	offset_sum <- rowSums(exp(offsetS))
	dataA <- data.frame(ysum=ysum, exp=design[,2], rep=design[,3], hic=design[,4], or2=design[,5], offset_sum=offset_sum)

	rownames(dataA) = sel.gene
	dataA_sel = dataA[intersect(which, sil.mutab[,1]),]

######### (strategy 1 : negative binomial) ### get all_betas from negative binomial model, needs optimization .....  
	if (method == "NB"){
		fit = glm.nb(ysum ~ exp + rep + hic + or2 + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100))
		
		###get mu_hat_glm
		mu_glm <- as.numeric(exp(design%*%as.matrix(fit$coefficients)))

		### NB: joint estimation of mu_hat and phi_hat	
		mu_hat <- get_muhat(y=y, offset=offsetS) 
		phi_hat_mle <- get_genewise_dispersion_mle_mu(y=y, mu=mu_hat, offset=offsetS, cores=cores)

		mu_hat_pre <- rep(0,length(mu_hat)) 
		while(sum(abs(mu_hat-mu_hat_pre)) > 10){
			###get phi, dispersion
			mu_hat_pre = mu_hat
			mu_hat <- get_mu_hat_mle(y=y, phi=phi_hat_mle, offset=offsetS, cores=cores)
			phi_hat_mle <- get_genewise_dispersion_mle_mu(y=y, mu=mu_hat, offset=offsetS, cores=cores)
			print(sum(abs(mu_hat-mu_hat_pre)))
		}

		###get sigma
		sigma <- getsigma(mu=mu_hat, mu_glm=mu_glm)
		
  		mu_post <- getmu_post_optimize(y=y, offset=offsetS, mu1=mu_hat, mu2=mu_glm, phi=phi_hat_mle, sigma = sigma, span=100, cores=cores)  		 
		mu_post_pre =  rep(0, length(mu_hat))
		while(sum(abs(mu_post-mu_post_pre)) > 10){
  			mu_post_pre <- mu_post
  			phi_hat_mle <- get_genewise_dispersion_mle_mu(y=y, mu=mu_post, offset=offsetS, cores=cores)
  			mu_post <- getmu_post_optimize(y=y, offset=offsetS, mu1=mu_hat, mu2=mu_glm, phi=phi_hat_mle, sigma = sigma, span=100, cores=cores)  
			print(sum(abs(mu_post-mu_post_pre)))
		}
		phi_hat = phi_hat_mle
	}


########## (strategy 2: poisson with log-linear) #### get all_betas from Poisson model
	if (method == "Poisson"){
		fit = glm(ysum ~ exp + rep + hic + or2 + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100),family=poisson)
	
		###get mu_hat_glm
		mu_glm <- as.numeric(exp(design%*%as.matrix(fit$coefficients)))
	
		### mu_hat for Poisson regression
		mu_hat <- rowSums(y)/rowSums(exp(offsetS))
	
		### get sigma
		sigma <- getsigma(mu=mu_hat, mu_glm=mu_glm)
	
		### get optimized mu 
		mu_post <- getmu_post_optimize(y=y, offset=offsetS, mu1=mu_hat, mu2=mu_glm, phi=0, sigma = sigma, span=100, cores=cores)
		phi_hat = rep(0, length(sel.gene))
	}
#############################################################
	
	### gene factor calculation ###	
	base = exp(fit$coefficients[1] + fit$coefficients[2] * as.numeric(gene.rep.expr[gid,3]) + fit$coefficients[3] * as.numeric(gene.rep.expr[gid,4]) + fit$coefficients[4] * as.numeric(gene.rep.expr[gid,5]) + fit$coefficients[5] * (as.numeric(as.character(gene.rep.expr[gid,6])) - 1))
        names(base) = gid 
	
	### assign correction factor based on base ############
	miss.gene = gid[gid %in% sel.gene==FALSE]
	or = order(base[sel.gene])
	sel.adj = data.table(data.frame(disp = phi_hat[or], p.base = (base[sel.gene])[or]))
	sel.adj$bin = ceiling(1:nrow(sel.adj)/50)
	disp.avg = sel.adj[, {dispAvg = median(disp); list(dispAvg=dispAvg)}, by='bin']
	setkey(disp.avg,bin)

	miss.data = base[miss.gene]
	miss.data[is.na(miss.data)] = 0
	range = findInterval(miss.data, sel.adj$p.base)
	range[range < 1] =1
	range[range > nrow(sel.adj)] = nrow(sel.adj)
	miss.disp = disp.avg[J(sel.adj[range,]$bin),]$dispAvg
	 
	### adjustment: dispersion, beta ###
	gene.phi = rep(0,length(gid))
	names(gene.phi) = gid
	gene.mu = gene.phi

	gene.phi[sel.gene] = phi_hat
	gene.phi[miss.gene] = miss.disp
	
	gene.mu[sel.gene] = mu_post

	offset = log(as.matrix(as.numeric(gene.prob[J(miss.gene,0),]$V1)) %*% t(as.matrix(as.numeric(mut.count[sample.name,2]))))
	offsetS <- offset - o_scale ### rescale 
	y = unname(count[miss.gene,])
	mu_post_miss <- getmu_post_optimize(y=y, offset=offsetS, mu1=rep(0,length(miss.gene)), mu2=base[miss.gene], phi=miss.disp, sigma = sigma, span=100, cores=cores)
	#gene.mu[miss.gene] = base[miss.gene]
	mu_post_miss[which(is.na(rowSums(offsetS)))] = base[miss.gene[which(is.na(rowSums(offsetS)))]] 
	gene.mu[miss.gene] = mu_post_miss

        rm(X)

        return(list( gene.mu, gene.prob, mut.count,gene.phi, o_scale))
}

# calculate background, and significance for SMG and TSG # 
calculate.bij<-function(gid, exome.gene, gene.mu, gene.prob, mut.count, gene.phi, sample.name, p,mut.context,nonsil.mutab, p_frameshift,o_scale)
{
	u = data.table(exome.gene)
	x = u[, sum(V1), by = 'gene']	
	gene.length = x$V1/3
	names(gene.length) = x$gene

	bij <- matrix(NA,nrow = length(gid), ncol = length(sample.name))
	rownames(bij) = gid
	colnames(bij) = sample.name
	stop.bij = bij

        res = rep(1,length(gid))
        names(res) = gid
        res.stop = res
        
	sampleVector = mut.count[sample.name,2]

  	res = rep(1,length(gid))
  	names(res) = gid
  	res.stop = res

 ##########
  	nosam <- aggregate(data=nonsil.mutab, Tumor_Sample_Barcode~ Ensembl_gene_id, function(x) length(unique(x))-1)
  	mut.matrix <- matrix(-1, nrow=length(gid), ncol =2)
  	rownames(mut.matrix) = gid
  	mut.matrix[nosam[,1],1] = nosam[,2]
  	nosam <- aggregate(data=nonsil.mutab[nonsil.mutab$Variant_Classification %in% c("Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del"),], Tumor_Sample_Barcode~ Ensembl_gene_id, function(x) length(unique(x))-1)
  	mut.matrix[nosam[,1],2] = nosam[,2]
 ############
	
	nonsil.p = as.matrix(gene.mu * ((gene.prob[J(gid,1),]$V1 + gene.length[gid]*p["indel"]) %*% t(as.matrix(sampleVector)))/ exp(o_scale))

	X = data.table(data.frame(gene = exome.gene[,1], p=p[as.numeric(exome.gene[,2])]* as.numeric(exome.gene[,4]), ind = as.numeric(exome.gene[,3])))
        gene.prob.2 = X[,sum(p), by = 'gene,ind']
        setkey(gene.prob.2,gene,ind)
	x = gene.prob.2[J(gid,2),]$V1
	x[is.na(x)] = 0
	stop.p = as.matrix(gene.mu * ((x + gene.length[gid]*p["indel"]*p_frameshift) %*% t(as.matrix(sampleVector))) / exp(o_scale))

	if(length(sample.name) < 2000){
        	med = "DFT-CF"
  	}else{
        	med = "RNA"
  	}

	for (i in 1:nrow(nonsil.p)){
		if ( i %% 100 == 0){cat(" ..", round(i/length(gid)*100),"% .. ",'\r'); flush.console()}
		gene = gid[i]
		if (gene.phi[gene] ==0){
			bij[i,] = exp(1)^(- nonsil.p[i,])	
		}else{	
			bij[i,] = (1/gamma(1))* (1/(1 + nonsil.p[i,] * gene.phi[gene]))^(1/gene.phi[gene])
		}
		if (mut.matrix[gene,1] >=0){
			res[gene] = 1-ppoibin(mut.matrix[gene,1], 1 - bij[i,],method = med)
		}
		if (mut.matrix[gene,2] >=0){
			if (gene.phi[gene]==0){
				stop.bg = exp(1)^(-stop.p[i,])
			}else{
				stop.bg = (1/gamma(1))* (1/(1 + stop.p[i,] * gene.phi[gene]))^(1/gene.phi[gene])
			}
			res.stop[gene] = 1-ppoibin(mut.matrix[gene,2], 1 - stop.bg,method = med)
		} 
	}
	
  	return(list(bij,res, res.stop, gene.length))
}


# calculate significance for OG activities # 
calculate.spatial <- function(gid, nonsilent.mutation.table, silent.mutation.table, exome, nosam, mut.context, p, p_inframe, gene.mu, mut.count, o_scale,og.num, bin.size, interval.size){
	### not consider in_frame mutation yet ###########        
	mutab = nonsilent.mutation.table[ (nonsilent.mutation.table[,"Variant_Classification"] %in% c("Missense_Mutation") & nonsilent.mutation.table[,"Variant_Type"] %in% c("Frame_shift","In_frame")==FALSE), ] #| nonsilent.mutation.table[,"Variant_Classification"] %in% c("In_Frame_Ins","In_Frame_Del"),]

	m = aggregate(Protein_position ~ Ensembl_gene_id, mutab, length)	
	use.gene = m[m[,2] >= og.num,1]

	mutab.filter= NULL
	mut.context.mod = data.table(data.frame( b = substr(mut.context[,3], 1, 3), a = substr(mut.context[,3], 4, 4), num = mut.context[,5]))
        setkey(mut.context.mod, b, a)

	
	l <- vector("list", 23) 
	
	for (i in 1:23){
		cat(" .. chromosome:",i ,"..",'\r'); flush.console()
                X = exome[[i]][as.character(exome[[i]][,7]) %in% use.gene,]
		if (nrow(X) ==0){next}
                an = as = p[mut.context.mod[J(X[,2],"1"),]$num]
                cn = cs = p[mut.context.mod[J(X[,2],"4"),]$num]
                gn = gs = p[mut.context.mod[J(X[,2],"3"),]$num]
                tn = ts = p[mut.context.mod[J(X[,2],"2"),]$num]
                an[X[,3] ==0 | X[,3] >=2] = 0
                cn[X[,4] ==0 | X[,4] >=2] = 0
                gn[X[,5] ==0 | X[,5] >=2] = 0
                tn[X[,6] ==0 | X[,6] >=2] = 0

                as[ X[,3] >=1 | is.na(as)] = 0
                cs[ X[,4] >=1 | is.na(cs)] = 0
                gs[ X[,5] >=1 | is.na(gs)] = 0
                ts[ X[,6] >=1 | is.na(ts)] = 0
		
		u = unlist(strsplit(as.character(X$V8),split = "/"))
                new.X = data.table(data.frame(nonsil.p = as.numeric(an + cn +gn +tn), sil.p = as.numeric(as + cs +gs + ts), gene = X[,7], pos = u[seq(1,length(u),2)], len = u[seq(2,length(u),2)]))

                new.X = new.X[, {v1 = sum(nonsil.p); v2 = sum(sil.p); list(v1=v1, v2=v2)}, by = 'gene,pos,len']
		l[[i]] = new.X
		X = data.table(exome[[i]])
		setkey(X,V1)
		sel.mutab = mutab[mutab[,2]==i,]
		mut = substr(sel.mutab$Context, 4,4)
		mut = as.numeric(chartr("1432", "3456", mut))
		sel.pos = data.frame(X[J(sel.mutab[,3]),])
		real.misense = which(sel.pos[col(sel.pos) == mut] == "1")
		sel.mutab[real.misense,]$Protein_position = as.character(sel.pos[real.misense,8])
		mutab.filter = rbind(mutab.filter, sel.mutab[real.misense,c("Chromosome","Protein_position","Ensembl_gene_id","Tumor_Sample_Barcode"),])
	}
	
	exome.pos = rbindlist(l)
	setkey(exome.pos,gene)


	# calculate genes	
	mutab = data.table(unique(mutab.filter))
	ac.count = mutab[, length(Protein_position), by = 'Ensembl_gene_id']

	sel.gene = ac.count[ac.count$V1 >=og.num,]$Ensembl_gene_id    
        mutab = data.frame(mutab[mutab$Ensembl_gene_id %in% sel.gene,])
	u = sapply(strsplit(as.character(mutab$Protein_position),split = "/"),unlist)
	mutab$Protein_position = u[1,]
	mutab$len = u[2,]
	
	########


	#g= "ENSG00000141510"
	#g="ENSG00000155657"
	cat ("\n",sprintf(paste0("Gene tested: ", length(sel.gene))),"\n")
	res = rep(1,length(gid))
        names(res) = gid
	res.pos = rep(NA, length(gid))
	names(res.pos) = gid


	for (i in 1:length(sel.gene)){
		if ( i %% 100 == 0){cat(" ..", round(i/length(sel.gene)*100),"% .. ",'\r'); flush.console()}
		g = sel.gene[i]
		m = mutab[mutab[,3]==g,]
		len = as.numeric(names(sort(table(m$len),decreasing=TRUE)[1]))	
        	ori.p = exome.pos[J(g),]
		s = as.numeric(as.character(ori.p$len)) == len & ori.p$pos != "-"
		ori.p = ori.p[s,]
		ori.p$pos = as.numeric(as.character(ori.p$pos))
		ori.p = ori.p[ori.p$pos <= as.numeric(as.character(ori.p$len)),]
		setkey(ori.p,pos)
		m = m[m$len == len,]

		og_correct = 1 - (gene.mu[g] * sum(ori.p$v1) * sum(mut.count[,2]) / exp(o_scale)) / nrow(m)
		if (og_correct < 0){og_correct = 0}
		
		ori.p$bin = ceiling(ori.p$pos/bin.size)
		p =rep(0, max(ori.p$bin))
		p[ori.p[, sum(v1), by= 'bin']$bin] = ori.p[, sum(v1), by= 'bin']$V1
		p = p/sum(p)
		
		n = rep(0, length(p))
		u = table(ori.p[J(as.numeric(m$Protein_position)),]$bin) * og_correct
		n[as.numeric(names(u))] = u
		na = which(round(n) >= 2)
		if (sum(n) < 1){next}	
	
		p.value = rep(1,length(p))
		for (j in na){
			p.value[j] = binom.test(round(n[j]), sum(round(n)), p = p[j], alternative = "greater")$p.value
		}

		
		#sig.pos = which(p.adjust(p.value, method = "fdr") < 0.05)
		sig.pos = na[which(p.adjust(p.value[na], method = "fdr") < 0.05)]
		if (length(sig.pos)==0){next}
		pos.new = sig.pos[n[sig.pos]>1]
		close.pos = which(diff(sig.pos) <=interval.size)
		
		for (j in close.pos){
			new.j = j
			group = NULL
			while(new.j %in% close.pos){
				group = c(group, sig.pos[new.j]:sig.pos[(new.j+1)])
				new.j = new.j+1
			}
			if (length(unique(m[m[,2] %in% group,4])) > 1){
				pos.new = c(pos.new,group)
			}
		}
		
		pos.new = unique(pos.new)
		res[g] = binom.test(sum(round(n[pos.new])), sum(round(n)), p = sum(p[pos.new]), alternative = "greater")$p.value
		if (length(pos.new) >0){
			pos.sort = sort(pos.new)
			dum.1 = c(0,diff(pos.sort))
                        dum.2 = which(dum.1 ==1)
                        if (length(dum.2) ==0){
                                dum.1[1:length(dum.1)] = ","
                        }else{  
                                dum.1[dum.2]="-"
                                dum.1[-dum.2] = ","             
                        }
                        dum.3 = paste(c(matrix(c(dum.1, pos.sort), 2, byrow = T)),sep="",collapse="")
                        a = gsub('-[^,]+-',"-",dum.3)
                        a = gsub('^,',"",a)

                        res.pos[g] = paste0("(",a,")/",len)
		}
	}
	return(cbind(res,res.pos))
}


# calculate pathway significance #
calculate.pathway <- function(bij, pathway.info, nonsil.mut.type,gene.symbols,gid){
        
	res = matrix(1,nrow=length(pathway.info),ncol=ncol(bij))  ### initial value
        res.multi = rep(1,length(pathway.info))   ###
        path.name = NULL

        if (ncol(bij) < 2000){
                med = "DFT-CF"
        }else{
                med = "RNA"
        }

	per.sample.cal <- function(x, mut, bij, sel.bij){
		mut.indicator = sum(mut[,2] == colnames(bij)[x]) -1
		sel.prob = sel.bij[,x]
		if (length(sel.prob) < 2000){
                	pval = 1-ppoibin(mut.indicator,1-sel.prob,method = "DFT-CF")
                }else{
                        pval = 1-ppoibin(mut.indicator,1-sel.prob,method = "RNA")
                }
		return(pval)
	}


        for (i in 1:length(pathway.info)){
                temp = unlist(strsplit(pathway.info[i],split="\t"))
                name = temp[1]
                path.name = c(path.name,name)

                sel.gene = intersect(rownames(gene.symbols)[match(temp[2:length(temp)],gene.symbols[,2])],gid)
                sel.bij = bij[sel.gene,]
                if (length(sel.bij) ==0){
                        next
                }
                mut = unique(nonsil.mut.type[nonsil.mut.type[,1] %in% sel.gene,c(1,8)])

                if (nrow(mut) ==0){
                        res[i,] = 1
                        res.multi[i] = 1
                        next
                }
                if(length(sel.bij) == ncol(bij)){
                        sel.bij = t(as.matrix(sel.bij))
                }
		
		
		res[i,] = unlist(mclapply(1:ncol(bij), per.sample.cal, mut=mut, bij=bij, sel.bij = sel.bij, mc.cores=1))


                sample.ind = sum(colnames(bij) %in% mut[,2]) - 1
                res.multi[i] = 1-ppoibin(sample.ind,1-exp(colSums(log(sel.bij),na.rm=T)),method = med)

        }
        rownames(res) = path.name
        colnames(res) = paste(colnames(bij),".p-value",sep="")
	pval.com = apply(res,1, function(x) combine.test(x, method = "fisher"))
    	res = cbind(res, pval.com, p.adjust(pval.com, method="fdr"))
	colnames(res)[(ncol(res)-1):ncol(res)] = c("MultiSample.p-value","MultiSample.q-value")

        names(res.multi) = path.name
        res.multi = cbind(res.multi,p.adjust(res.multi,method="fdr"))
        colnames(res.multi) = c("CrossSample.p-value","CrossSample.q-value")
        return(cbind(res,res.multi))
}







############################ sample - specific functions ######################
# get background count for silent / non-silent mutations #
preprocess.BM.sample <-function(X, gene.nonsil, gene, mut.context, id.to.name)
{
	gene = id.to.name[gene,1]
	gene.nonsil = id.to.name[gene.nonsil,1]

        res = X[, sum(V5), by = 'V3,V4,V2']
	rm(X)
	res[, sil := ifelse(V2 %in% gene, 1, 0)]
	res[, nonsil := ifelse(V2 %in% gene.nonsil, 1, 0)]
	
        sil = res[res$V4 ==0,]
        nonsil = res[res$V4 > 0,]


	###############################
	sil.type.const = rbind(rep(0,96,),rep(0,96))
        colnames(sil.type.const) = 1:96
	nonsil.type.const = sil.type.const

	# for silent 
	sil.summary = subset(sil[,sum(V1),by = 'V3,sil'],sil==1)	
	nonsil.summary = subset(nonsil[,sum(V1),by = 'V3,sil'], sil==1)

        e = sil.summary[match(colnames(sil.type.const),sil.summary$V3),]$V1
        f = nonsil.summary[match(colnames(sil.type.const),nonsil.summary$V3),]$V1
        e[is.na(e)] = 0
        f[is.na(f)] = 0

        sil.type.const = sil.type.const + rbind(e,f)

	# for non-silent
	sil.summary = subset(sil[,sum(V1),by = 'V3,nonsil'], nonsil==1)
        nonsil.summary = subset(nonsil[,sum(V1),by = 'V3,nonsil'],nonsil ==1)

        e = sil.summary[match(colnames(nonsil.type.const),sil.summary$V3),]$V1
        f = nonsil.summary[match(colnames(nonsil.type.const),nonsil.summary$V3),]$V1
        e[is.na(e)] = 0
        f[is.na(f)] = 0

        nonsil.type.const = nonsil.type.const + rbind(e,f)
	
        return(list(nonsil.type.const,sil.type.const))
}

### read in sample exomes ##
exome.sample.load <- function(beds)
{
	all.files <- paste(beds,list.files(path = beds,pattern = ".txt"),sep="/")
	mylist <- lapply(all.files, fread)
	exome.gene.sample <- rbindlist(mylist)
	#exome.gene.sample <- fread(beds)
	rm(mylist)
	return(exome.gene.sample)
}


# fit background model with negative-binomial regression #
fit.model.sample <- function(X,nonsil.mutab,sil.mutab,sample.name,nosam,gene.rep.expr,gid,p,S, which, method, id.to.name, cores){
	p.reshape = data.table(p = p[1:96], id = as.integer(names(p)[1:96]), key="id")
	X[,V3 := p.reshape[J(X$V3),p] * V5]

        gene.prob = X[,{prob = sum(V3); list(prob=prob)}, by = 'V2,V4,V1']
        setkey(gene.prob,V2,V4,V1)
 

        # sample count # 
	mut.count = aggregate(Start_Position ~ Tumor_Sample_Barcode, rbind(sil.mutab,nonsil.mutab),length)
	rownames(mut.count) = mut.count[,1]
	sample.length = X[,{len=sum(V5);list(len=len)}, by = 'V1']
	setkey(sample.length,V1)
	mut.count[,2] = mut.count[,2]/(sample.length[J(id.to.name[mut.count[,1],1]),]$len/3)
	
	rm(X)

        # for all genes with silent mutations 
        sel.gene = intersect(rownames(gene.rep.expr), sil.mutab[,1])
	offset = matrix(0, nrow=length(sel.gene),ncol=length(sample.name))
	
	for (i in 1:length(sample.name)){	
		x = gene.prob[J(id.to.name[sel.gene,1],0,id.to.name[sample.name[i],1]),]$prob
		x[is.na(x)]=0
		x = x * mut.count[sample.name[i],2]
		offset[,i] = x
	}  ##### offset in normal scale ##########


        count = matrix(0,nrow = length(gid), ncol = S)
        rownames(count) = gid
        colnames(count) = sample.name
        sil.count = aggregate(Start_Position ~ Ensembl_gene_id + Tumor_Sample_Barcode, sil.mutab,length)
        for (i in 1:nrow(sil.count)){
                count[as.character(sil.count[i,1]),as.character(sil.count[i,2])] = sil.count[i,3]
        }
        y = unname(count[sel.gene,])   #### y 
        m= data.frame(exp = as.numeric(gene.rep.expr[sel.gene,3]), rep = as.numeric(gene.rep.expr[sel.gene,4]), hic = as.numeric(gene.rep.expr[sel.gene,5]), or = factor(gene.rep.expr[sel.gene,6], levels =1:2))
        design = model.matrix(~ exp + rep + hic + or, m) ### design
	
	ysum <- apply(y, 1, sum)

	o_scale <- max(offset)
        offset <- offset/o_scale ### rescale 
	
	offset_sum <- rowSums(offset)
	dataA <- data.frame(ysum=ysum, exp=design[,2], rep=design[,3], hic=design[,4], or2=design[,5], offset_sum=offset_sum)

	rownames(dataA) = sel.gene
	dataA_sel = dataA[intersect(which, sil.mutab[,1]),]


######### (strategy 1 : negative binomial) ### get all_betas from negative binomial model, needs optimization .....  
	if (method == "NB"){
		fit = glm.nb(ysum ~ exp + rep + hic + or2 + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100))

		###get mu_hat_glm
		mu_glm <- as.numeric(exp(design%*%as.matrix(fit$coefficients)))

		### NB: joint estimation of mu_hat and phi_hat	
		mu_hat <- get_muhat_sample(y=y, offset=offset) 
		phi_hat_mle <- get_genewise_dispersion_mle_mu_sample(y=y, mu=mu_hat, offset=offset, cores=cores)

		mu_hat_pre <- rep(0,length(mu_hat)) 
		while(sum(abs(mu_hat-mu_hat_pre)) > 10){
			###get phi, dispersion
			mu_hat_pre = mu_hat
			mu_hat <- get_mu_hat_mle_sample(y=y, phi=phi_hat_mle, offset=offset, cores=cores)
			phi_hat_mle <- get_genewise_dispersion_mle_mu_sample(y=y, mu=mu_hat, offset=offset,cores=cores)
			print(sum(abs(mu_hat-mu_hat_pre)))
		}

		###get sigma
		sigma <- getsigma(mu=mu_hat, mu_glm=mu_glm)
		
  		mu_post <- getmu_post_optimize_sample(y=y, offset=offset, mu1=mu_hat, mu2=mu_glm, phi=phi_hat_mle, sigma = sigma, span=100, cores=cores)  		 
		mu_post_pre =  rep(0, length(mu_hat))
		while(sum(abs(mu_post-mu_post_pre)) > 10){
  			mu_post_pre <- mu_post
  			phi_hat_mle <- get_genewise_dispersion_mle_mu_sample(y=y, mu=mu_post, offset=offset, cores=cores)
  			mu_post <- getmu_post_optimize_sample(y=y, offset=offset, mu1=mu_hat, mu2=mu_glm, phi=phi_hat_mle, sigma = sigma, span=100, cores=cores)  
			print(sum(abs(mu_post-mu_post_pre)))
		}
		phi_hat = phi_hat_mle
	}


########## (strategy 2: poisson with log-linear) #### get all_betas from Poisson model
	if (method == "Poisson"){
		fit = glm(ysum ~ exp + rep + hic + or2 + offset(log(offset_sum)), data = dataA_sel,control=glm.control(maxit=100),family=poisson)
	
		###get mu_hat_glm
		mu_glm <- as.numeric(exp(design%*%as.matrix(fit$coefficients)))
	
		### mu_hat for Poisson regression
		mu_hat <- rowSums(y)/rowSums(offset)
	
		### get sigma
		sigma <- getsigma(mu=mu_hat, mu_glm=mu_glm)
	
		### get optimized mu 
		mu_post <- getmu_post_optimize_sample(y=y, offset=offset, mu1=mu_hat, mu2=mu_glm, phi=0, sigma = sigma, span=100,cores=cores)
		phi_hat = rep(0, length(sel.gene))
	}
#############################################################
	
	### gene factor calculation ###	
	base = exp(fit$coefficients[1] + fit$coefficients[2] * as.numeric(gene.rep.expr[gid,3]) + fit$coefficients[3] * as.numeric(gene.rep.expr[gid,4]) + fit$coefficients[4] * as.numeric(gene.rep.expr[gid,5]) + fit$coefficients[5] * (as.numeric(as.character(gene.rep.expr[gid,6])) - 1))
        names(base) = gid 
	
	### assign correction factor based on base ############
	miss.gene = gid[gid %in% sel.gene==FALSE]
	or = order(base[sel.gene])
	sel.adj = data.table(data.frame(disp = phi_hat[or], p.base = (base[sel.gene])[or]))
	sel.adj$bin = ceiling(1:nrow(sel.adj)/50)
	disp.avg = sel.adj[, {dispAvg = median(disp); list(dispAvg=dispAvg)}, by='bin']
	setkey(disp.avg,bin)

	miss.data = base[miss.gene]
	miss.data[is.na(miss.data)] = 0
	range = findInterval(miss.data, sel.adj$p.base)
	range[range < 1] =1
	range[range > nrow(sel.adj)] = nrow(sel.adj)
	miss.disp = disp.avg[J(sel.adj[range,]$bin),]$dispAvg
	 
	### adjustment: dispersion, beta ###
	gene.phi = rep(0,length(gid))
	names(gene.phi) = gid
	gene.mu = gene.phi

	gene.phi[sel.gene] = phi_hat
	gene.phi[miss.gene] = miss.disp

	gene.mu[sel.gene] = mu_post
	

	offset = matrix(0, nrow=length(miss.gene),ncol=length(sample.name))
        for (i in 1:length(sample.name)){
                x = gene.prob[J(id.to.name[miss.gene,1],0,id.to.name[sample.name[i],1]),]$prob
                x[is.na(x)]=0
                x = x * mut.count[sample.name[i],2]
                offset[,i] = x
        }  ##### offset in normal scale ##########

        offset <- offset/o_scale 
        y = unname(count[miss.gene,])
        mu_post_miss <- getmu_post_optimize_sample(y=y, offset=offset, mu1=rep(0,length(miss.gene)), mu2=base[miss.gene], phi=miss.disp, sigma = sigma, span=100, cores=cores)
        gene.mu[miss.gene] = mu_post_miss

	#gene.mu[miss.gene] = base[miss.gene]


        return(list( gene.mu, gene.prob, mut.count,gene.phi, sample.length,log(o_scale)))
}

# calculate background, and significance for SMG and TSG # 
calculate.bij.sample<-function(gid, X, gene.mu, gene.prob, mut.count, gene.phi, sample.name, p,mut.context,nonsil.mutab, p_frameshift, sample.length, id.to.name, o_scale)
{
	gene.length = X[, {len=sum(V5)/3; list(len=len)}, by = 'V1,V2']	
	setkey(gene.length, V1, V2) 

	bij <- matrix(NA,nrow = length(gid), ncol = length(sample.name))
	rownames(bij) = gid
	colnames(bij) = sample.name
	stop.bij = bij

        res = rep(1,length(gid))
        names(res) = gid
        res.stop = res
        
	sampleVector = mut.count[sample.name,2]

  	res = rep(1,length(gid))
  	names(res) = gid
  	res.stop = res

 ##########
  	nosam <- aggregate(data=nonsil.mutab, Tumor_Sample_Barcode~ Ensembl_gene_id, function(x) length(unique(x))-1)
  	mut.matrix <- matrix(-1, nrow=length(gid), ncol =2)
  	rownames(mut.matrix) = gid
  	mut.matrix[nosam[,1],1] = nosam[,2]
  	nosam <- aggregate(data=nonsil.mutab[nonsil.mutab$Variant_Classification %in% c("Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del"),], Tumor_Sample_Barcode~ Ensembl_gene_id, function(x) length(unique(x))-1)
  	mut.matrix[nosam[,1],2] = nosam[,2]
 ############

	gene.prob.nonsil = gene.prob[gene.prob$V4 >0,]
	gene.nonsil.sum = gene.prob[,{p=sum(prob); list(p=p)}, by ='V1,V2']
	setkey(gene.nonsil.sum, V1,V2)
	setkey(gene.prob, V1,V2,V4)

	nonsil.mean = stop.mean = matrix(0, nrow=length(gid), ncol= length(sample.name))
	for (i in 1:length(sample.name)){
		x = gene.mu * (gene.nonsil.sum[J(id.to.name[sample.name[i],1], id.to.name[gid,1]),]$p + gene.length[J(id.to.name[sample.name[i],1], id.to.name[gid,1]),]$len * p["indel"]) * sampleVector[i] / exp(o_scale)
		x[is.na(x)] = 0
		nonsil.mean[,i] = x
				
		x = gene.mu * (gene.prob[J(id.to.name[sample.name[i],1], id.to.name[gid,1],2),]$prob + gene.length[J(id.to.name[sample.name[i],1], id.to.name[gid,1]),]$len * p["indel"] * p_frameshift) * sampleVector[i] / exp(o_scale)
		x[is.na(x)] = 0
		stop.mean[,i] = x	
	}


	if(length(sample.name) < 2000){
        	med = "DFT-CF"
  	}else{
        	med = "RNA"
  	}

	for (i in 1:nrow(nonsil.mean)){
		if ( i %% 100 == 0){cat(" ..", round(i/length(gid)*100),"% .. ",'\r'); flush.console()}
		gene = gid[i]
		if (gene.phi[gene] ==0){
			bij[i,] = exp(1)^(- nonsil.mean[i,])	
		}else{	
			bij[i,] = (1/gamma(1))* (1/(1 + nonsil.mean[i,] * gene.phi[gene]))^(1/gene.phi[gene])
		}
		if (mut.matrix[gene,1] >=0){
			res[gene] = 1-ppoibin(mut.matrix[gene,1], 1 - bij[i,],method = med)
		}
		if (mut.matrix[gene,2] >=0){
			if (gene.phi[gene]==0){
				stop.bg = exp(1)^(-stop.mean[i,])
			}else{
				stop.bg = (1/gamma(1))* (1/(1 + stop.mean[i,] * gene.phi[gene]))^(1/gene.phi[gene])
			}
			res.stop[gene] = 1-ppoibin(mut.matrix[gene,2], 1 - stop.bg,method = med)
		} 
	}
	
  	return(list(bij,res, res.stop, gene.length))
}



## calculate significance for OG activities # 
calculate.spatial.sample <- function(gid, nonsilent.mutation.table, silent.mutation.table, exome, nosam, mut.context, p, p_inframe, gene.mu, mut.count, o_scale,og.num, bin.size, interval.size){
	### not consider in_frame mutation yet ###########        
	mutab = nonsilent.mutation.table[ (nonsilent.mutation.table[,"Variant_Classification"] %in% c("Missense_Mutation") & nonsilent.mutation.table[,"Variant_Type"] %in% c("Frame_shift","In_frame")==FALSE), ] #| nonsilent.mutation.table[,"Variant_Classification"] %in% c("In_Frame_Ins","In_Frame_Del"),]

	m = aggregate(Protein_position ~ Ensembl_gene_id, mutab, length)	
	use.gene = m[m[,2] >= og.num,1]

	mutab.filter= NULL
	mut.context.mod = data.table(data.frame( b = substr(mut.context[,3], 1, 3), a = substr(mut.context[,3], 4, 4), num = mut.context[,5]))
        setkey(mut.context.mod, b, a)


	l <- vector("list", 23) 
 
	
	for (i in 1:23){
		cat(" .. chromosome:",i ,"..",'\r'); flush.console()
                X = exome[[i]][as.character(exome[[i]][,7]) %in% use.gene,]
		if (nrow(X) ==0){next}
                an = as = p[mut.context.mod[J(X[,2],"1"),]$num]
                cn = cs = p[mut.context.mod[J(X[,2],"4"),]$num]
                gn = gs = p[mut.context.mod[J(X[,2],"3"),]$num]
                tn = ts = p[mut.context.mod[J(X[,2],"2"),]$num]
                an[X[,3] ==0 | X[,3] >=2] = 0
                cn[X[,4] ==0 | X[,4] >=2] = 0
                gn[X[,5] ==0 | X[,5] >=2] = 0
                tn[X[,6] ==0 | X[,6] >=2] = 0

                as[ X[,3] >=1 | is.na(as)] = 0
                cs[ X[,4] >=1 | is.na(cs)] = 0
                gs[ X[,5] >=1 | is.na(gs)] = 0
                ts[ X[,6] >=1 | is.na(ts)] = 0
		
		u = unlist(strsplit(as.character(X$V8),split = "/"))
                new.X = data.table(data.frame(nonsil.p = as.numeric(an + cn +gn +tn), sil.p = as.numeric(as + cs +gs + ts), gene = X[,7], pos = u[seq(1,length(u),2)], len = u[seq(2,length(u),2)]))

                new.X = new.X[, {v1 = sum(nonsil.p); list(v1=v1)}, by = 'gene,pos,len']
		l[[i]] = new.X
		X = data.table(exome[[i]])
		setkey(X,V1)
		sel.mutab = mutab[mutab[,2]==i,]
		mut = substr(sel.mutab$Context, 4,4)
		mut = as.numeric(chartr("1432", "3456", mut))
		sel.pos = data.frame(X[J(sel.mutab[,3]),])
		real.misense = which(sel.pos[col(sel.pos) == mut] == "1")
		sel.mutab[real.misense,]$Protein_position = as.character(sel.pos[real.misense,8])
		mutab.filter = rbind(mutab.filter, sel.mutab[real.misense,c("Chromosome","Protein_position","Ensembl_gene_id","Tumor_Sample_Barcode"),])
	}

	exome.pos = rbindlist(l)
	setkey(exome.pos,gene)


	# calculate genes	
	mutab = data.table(unique(mutab.filter))
	ac.count = mutab[, length(Protein_position), by = 'Ensembl_gene_id']

	sel.gene = ac.count[ac.count$V1 >= og.num,]$Ensembl_gene_id    
        mutab = data.frame(mutab[mutab$Ensembl_gene_id %in% sel.gene,])
	u = sapply(strsplit(as.character(mutab$Protein_position),split = "/"),unlist)
	mutab$Protein_position = u[1,]
	mutab$len = u[2,]
	
	########

	cat ("\n",sprintf(paste0("Gene tested: ", length(sel.gene))),"\n")
	res = rep(1,length(gid))
        names(res) = gid
	res.pos = rep(NA, length(gid))
	names(res.pos) = gid


	for (i in 1:length(sel.gene)){
		if ( i %% 100 == 0){cat(" ..", round(i/length(sel.gene)*100),"% .. ",'\r'); flush.console()}
		g = sel.gene[i]
		m = mutab[mutab[,3]==g,]
		len = as.numeric(names(sort(table(m$len),decreasing=TRUE)[1]))	
        	ori.p = exome.pos[J(g),]
		s = as.numeric(as.character(ori.p$len)) == len & ori.p$pos != "-"
		ori.p = ori.p[s,]
		ori.p$pos = as.numeric(as.character(ori.p$pos))
		ori.p = ori.p[ori.p$pos <= as.numeric(as.character(ori.p$len)),]
		setkey(ori.p,pos)
		m = m[m$len == len,]

		og_correct = 1 - (gene.mu[g] * sum(ori.p$v1) * sum(mut.count[,2]) / exp(o_scale)) / nrow(m)
		if (og_correct < 0){og_correct = 0}
		
		ori.p$bin = ceiling(ori.p$pos/bin.size)
		p =rep(0, max(ori.p$bin))
		p[ori.p[, sum(v1), by= 'bin']$bin] = ori.p[, sum(v1), by= 'bin']$V1
		p = p/sum(p)
		
		n = rep(0, length(p))
		u = table(ori.p[J(as.numeric(m$Protein_position)),]$bin) * og_correct
		n[as.numeric(names(u))] = u
		na = which(round(n) >= 2)   ######### at least 2 mutations at this position
		if (sum(n) < 1){next}	
	
		p.value = rep(1,length(p))
		for (j in na){
			p.value[j] = binom.test(round(n[j]), sum(round(n)), p = p[j], alternative = "greater")$p.value
		}

		
		sig.pos = na[which(p.adjust(p.value[na], method = "fdr") < 0.05)]
		if (length(sig.pos)==0){next}
		pos.new = sig.pos[n[sig.pos]>1]
		close.pos = which(diff(sig.pos) <=interval.size)
		
		for (j in close.pos){
			new.j = j
			group = NULL
			while(new.j %in% close.pos){
				group = c(group, sig.pos[new.j]:sig.pos[(new.j+1)])
				new.j = new.j+1
			}
			if (length(unique(m[m[,2] %in% group,4])) > 1){
				pos.new = c(pos.new,group)
			}
		}
		
		pos.new = unique(pos.new)
		res[g] = binom.test(sum(round(n[pos.new])), sum(round(n)), p = sum(p[pos.new]), alternative = "greater")$p.value
		if (length(pos.new) >0){
			pos.sort = sort(pos.new)
			dum.1 = c(0,diff(pos.sort))
			dum.2 = which(dum.1 ==1)
			if (length(dum.2) ==0){
				dum.1[1:length(dum.1)] = ","
			}else{
				dum.1[dum.2]="-"
				dum.1[-dum.2] = ","		
			}
			dum.3 = paste(c(matrix(c(dum.1, pos.sort), 2, byrow = T)),sep="",collapse="")
			a = gsub('-[^,]+-',"-",dum.3)
			a = gsub('^,',"",a)

			res.pos[g] = paste0("(",a,")/",len)
		}
	}
	return(cbind(res,res.pos))
}







################## main function ##
driver.detection <- function(maf.file, exome.file, exome.gene, covar.file, mutation.context, fraction, subtype.file, pathway.file, cgc, beds, model, cores,og.num, bin.size, interval.size) {
  	
	# load data
	mut.context = read.table(mutation.context)
	
	gene.rep.expr <- read.table(covar.file,header=T,row.names=2,sep="\t")

  	exome.ptr <- load(exome.file)
	exome <- get(exome.ptr)
  	rm(exome.ptr)
	
	exome.ptr <- load(exome.gene)
	exome.gene <-  get(exome.ptr)
  	exome.gene$tag = mut.context[match(exome.gene$tag,mut.context[,3]),5]

  	gid <- intersect(rownames(gene.rep.expr)[gene.rep.expr[,1] %in% c(1:22,"X")],exome.gene$gene)

	rm(exome.ptr)
	
	# read in maf file
	maf.table = read.table(maf.file,header=T,sep="\t",stringsAsFactors=FALSE)
	maf.table = maf.table[maf.table[,1] %in% gid,]
	
	cat(sprintf(paste0("Total gene analyzed: ",length(gid))),"\n")

	nonsilent.mutation.table = maf.table[maf.table$Sil==1,]
	silent.mutation.table = maf.table[maf.table$Sil==2,]

	nonsilent.mutation.table = get.mutation.exome(exome,nonsilent.mutation.table)
	silent.mutation.table = get.mutation.exome(exome,silent.mutation.table) 
 
	nosam <- aggregate(data=nonsilent.mutation.table, Tumor_Sample_Barcode~ Ensembl_gene_id, function(x) length(unique(x)))	
	mut.gene = nosam[with(nosam,order(Tumor_Sample_Barcode)),1]
	nonsilent_passenger_gene=gid[gid %in% mut.gene[round(length(mut.gene)* fraction):length(mut.gene)] ==FALSE]

	x=nonsilent.mutation.table[ nonsilent.mutation.table[,1] %in% nonsilent_passenger_gene, 4 ]
        p_inframe=sum(x=="In_frame")/(sum(x=="In_frame")+ sum(x=="Frame_shift") )
        p_frameshift=sum(x=="Frame_shift")/(sum(x=="In_frame")+ sum(x=="Frame_shift") )
        if(p_inframe==0 | p_frameshift==0 | sum(x=="In_frame")+ sum(x=="Frame_shift") <=5 ){
                p_inframe= 1/3
                p_frameshift= 2/3
        }
	
	nonsil.mutab = nonsilent.mutation.table
	sil.mutab = silent.mutation.table

	sample.name=unique(c(as.character(sil.mutab[,8]),as.character(nonsil.mutab[,8])))     ### sample id
	S = length(sample.name)


 	 #####  calculate p_{i}, i=1,2...,96 and selection bias r ##############################################################################  
	nonsil.mut.type.sampl.sum <- mut.type.converter(nonsilent_passenger_gene,nonsil.mutab, mut.context)
  	sil.mut.type.sampl.sum <- mut.type.converter(gid,sil.mutab, mut.context)
	
	# sample info 
	if (beds != "NA"){
		exome.gene.sample <- exome.sample.load(beds)
		id.to.name = read.table(paste(beds,"/sample_gene_id",sep=""),row.names=1)
		res=preprocess.BM.sample(exome.gene.sample, nonsilent_passenger_gene, gid, mut.context,id.to.name)
		nonsil.type.const=res[[1]]
		sil.type.const = res[[2]]
		gid = intersect(rev(rownames(id.to.name))[match(unique(exome.gene.sample$V2),id.to.name[rev(rownames(id.to.name)),1])],gid)
	}else{
  		nonsil.type.const=preprocess.BM(exome.gene, nonsilent_passenger_gene, mut.context)  
	  	sil.type.const=preprocess.BM(exome.gene, gid,mut.context) 
	}

  	temp =  fit.background.model(nonsil.mutab[nonsil.mutab[,1] %in% nonsilent_passenger_gene,], nonsil.mut.type.sampl.sum, sil.mut.type.sampl.sum, nonsil.type.const, sil.type.const,mut.context)
	p = temp[[1]]

	##### fit possion regression model #####

	gene.rep.expr = as.matrix(gene.rep.expr,stringsAsFactors=FALSE)

	which = rownames(gene.rep.expr)[!is.na(gene.rep.expr[,3]) & !is.na(gene.rep.expr[,4]) & !is.na(gene.rep.expr[,5])]

	gene.rep.expr[is.na(gene.rep.expr[,3]),3] = mean(as.numeric(gene.rep.expr[,3]),na.rm=T)
	gene.rep.expr[is.na(gene.rep.expr[,4]),4] = mean(as.numeric(gene.rep.expr[,4]),na.rm=T)
	gene.rep.expr[is.na(gene.rep.expr[,5]),5] = mean(as.numeric(gene.rep.expr[,5]),na.rm=T)	 
		
	m = prcomp(apply(gene.rep.expr[,3:5],2,as.numeric),scale=T)
	gene.rep.expr[,3:5] = m$x

	cat(sprintf(paste("Start background mutation rate estimation"," ...", sep="")),"\n")

### sample info
	if (beds != "NA"){
		res = fit.model.sample(exome.gene.sample,nonsil.mutab,sil.mutab,sample.name,nosam, gene.rep.expr, gid, p, S, which,model, id.to.name, cores)	
		gene.mu = res[[1]]
		gene.prob = res[[2]]
		mut.count = res[[3]]
		gene.phi = res[[4]]
		sample.length = res[[5]]	
		o_scale = res[[6]]

		bg <- calculate.bij.sample(gid, exome.gene.sample, gene.mu, gene.prob, mut.count, gene.phi, sample.name, p,mut.context, nonsil.mutab, p_frameshift, sample.length, id.to.name,o_scale)
		bij = bg[[1]]
		res = bg[[2]] 
		tsg.p = bg[[3]]
		smg.p = cbind(as.character(gene.rep.expr[match(names(res),rownames(gene.rep.expr)),2]),res)	
		gene.length = bg[[4]]
	
	}else{
		res = fit.model(exome.gene,nonsil.mutab,sil.mutab,sample.name,nosam, gene.rep.expr, gid, p, S, which, model,cores)	

		gene.mu = res[[1]]
		gene.prob = res[[2]]
		mut.count = res[[3]]
		gene.phi = res[[4]]
		o_scale = res[[5]]

  		##### calculate bij  #####
  	 	bg <- calculate.bij(gid, exome.gene, gene.mu, gene.prob, mut.count, gene.phi, sample.name, p,mut.context, nonsil.mutab, p_frameshift, o_scale)	
		bij = bg[[1]]
		res = bg[[2]] 
		tsg.p = bg[[3]]
		smg.p = cbind(as.character(gene.rep.expr[match(names(res),rownames(gene.rep.expr)),2]),res)	
		gene.length = bg[[4]]
	}

	#### p-value from Recurrence ####
	cat(sprintf(paste("Gene-specific mutation rate calculation is done"," ...",sep="")),"\n")
	
	#tag = unlist(strsplit(rev(unlist(strsplit(maf.file, "/")))[1],"-clean2."))[1]
	#save(bij,file=paste("res/",tag,".Rdata",sep=""))

	#### p from spatial calculation, classify as oncogene
	cat(sprintf (paste("Start spatial calculation"," ...",sep="")),"\n")
	
	if (beds != "NA"){
		og.p <- calculate.spatial.sample(gid, nonsilent.mutation.table, silent.mutation.table, exome, nosam, mut.context, p, p_inframe, gene.mu, mut.count, o_scale,og.num,bin.size, interval.size)   ##### To Do 
	}else{
		og.p <- calculate.spatial(gid, nonsilent.mutation.table, silent.mutation.table, exome, nosam, mut.context, p, p_inframe, gene.mu, mut.count, o_scale,og.num, bin.size, interval.size) 
	}


	#### Pathway analysis ####
	cat(sprintf (paste("Driver gene calculation is done","; Start pathway analysis ...",sep="")),"\n")
	pathway.info = readLines(pathway.file)		
	path.pp <- calculate.pathway(bij,pathway.info,nonsilent.mutation.table,gene.rep.expr,gid)
	
	#### output ###
	mutab = data.table(rbind(nonsilent.mutation.table,silent.mutation.table))
	gene.count = mutab[,length(Start_Position),by = 'Ensembl_gene_id,Sil']
	sample.count = mutab[,length(unique(Tumor_Sample_Barcode)),by = 'Ensembl_gene_id,Sil']
	setkey(gene.count, Ensembl_gene_id,Sil)
	setkey(sample.count,Ensembl_gene_id,Sil)
	
	count = cbind(gene.count[J(gid,2),]$V1, sample.count[J(gid,2),]$V1, gene.count[J(gid,1),]$V1, sample.count[J(gid,1),]$V1)
	count[is.na(count)] = 0

	info = rep("",nrow(smg.p))
	info[smg.p[,1] %in% cgc[,1]] = "CGC"
	gene.pp = cbind(rownames(smg.p),smg.p[,1],info,count, smg.p[,2],tsg.p, og.p)

	#gene.pp = cbind(gene.pp,p.adjust(as.numeric(gene.pp[,8]),method="fdr"), p.adjust(as.numeric(gene.pp[,9]),method="fdr"),p.adjust(as.numeric(gene.pp[,10]),method="fdr"),apply(gene.pp, 1, function(x) combine.test(c(as.numeric(x[8]),as.numeric(x[9]),as.numeric(x[10])),method="fisher")))
	gene.pp = cbind(gene.pp,p.adjust(as.numeric(gene.pp[,8]),method="bonferroni"), p.adjust(as.numeric(gene.pp[,9]),method="bonferroni"),p.adjust(as.numeric(gene.pp[,10]),method="bonferroni"),apply(gene.pp, 1, function(x) combine.test(c(as.numeric(x[8]),as.numeric(x[9]),as.numeric(x[10])),method="fisher")))
	
	gene.tag = matrix(" ",nrow=nrow(gene.pp),ncol=3)
	gene.tag[as.numeric(gene.pp[,12])<0.05,1] = "SMG"
	gene.tag[as.numeric(gene.pp[,13])<0.05,2] = "TSG"
	gene.tag[as.numeric(gene.pp[,14])<0.05,3] = "OG"
	q.tag = apply(gene.tag, 1, function(x) paste(x,collapse= ","))
	q.tag = gsub(" ,","",q.tag)
	q.tag = gsub(" $","",q.tag)
	q.tag = gsub(",$","",q.tag)
	
	gene.pp = cbind(gene.pp, q.tag)


	colnames(gene.pp) = c("Ensembl_gene_id","Gene_symbol","CGC","Mut_silent","Sample_silent","Mut_nonsilent","Sample_nonsilent","P-value_SMG", "P-value_TSG","P-value_OG","OG_clusters_aa_pos","Q-value_SMG","Q-value_TSG","Q-value_OG", "Fisher_combined_P","Tag(Q<0.05)")

	gene.pp = gene.pp[order(as.numeric(gene.pp[,15])),]

  	rm(gene.rep.expr)
	rm(nonsilent.mutation.table)
	
	return(list(gene.pp, path.pp))
}


################################
####### run & write output #####
###############################

tag = unlist(strsplit(rev(unlist(strsplit(maf.file, "/")))[1],"-clean2."))[1]
file.gene = paste(output,"/",tag,"_",method,"_sig-genes.tsv",sep="")
file.pathway = paste(output,"/",tag,"_",method,"_sig-pathways.tsv",sep="")
cat(sprintf (paste("Start analysis on ", tag, sep="")),"\n")

res = driver.detection(maf.file, exome.file, exome.gene, covar.file, mutation.context, fraction, subtype.file, pathway.file, cgc, beds, method, cores, og.num, bin.size, interval.size)
write.table(res[[1]], file = file.gene, row.names=F,col.names=T,sep="\t",quote=F)
write.table(res[[2]], file = file.pathway, row.names=T, col.names=T, sep="\t",quote=F)
cat(sprintf (paste("Outputs are: ", file.gene, "; ", file.pathway, sep="")),"\n")

