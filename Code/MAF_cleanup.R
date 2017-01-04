args <- commandArgs(trailingOnly=TRUE)

maf.file <- args[1]
output <- args[2]



maf.cleanup <- function(maf.file,output){

        maf.table <- read.delim(maf.file, header=TRUE, sep="\t", stringsAsFactors=FALSE, skip=1)
        maf.table = as.matrix(apply(maf.table,2,as.character))
        maf.table <- maf.table[!is.na(maf.table[,5]),] #### remove positions without chromosome info.

        type.col <- which(colnames(maf.table)=="Variant_Type")
        class.col <- which(colnames(maf.table)=="Variant_Classification")
        nonsilent.maf.table = maf.table[maf.table[,class.col] %in% c("Missense_Mutation","Nonsense_Mutation","Splice_Site","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Ins","In_Frame_Del","Nonstop_Mutation","Translation_Start_Site"),]

        silent.maf.table <- maf.table[maf.table[,class.col] %in% c("Silent") & (maf.table[,type.col] %in% c("INS", "DEL") ==FALSE) ,]   ### table of nonsilent mutations

########### clean -up  MAF file #######
        which.indels <- which(nonsilent.maf.table[,type.col] %in% c("INS", "DEL"))
        if (length(which.indels)>0) {
                indel.labels <- nonsilent.maf.table[which.indels, class.col]
                if (sum(indel.labels %in% c("Frame_Shift_Del", "Frame_Shift_Ins","Nonsense_Mutation","Splice_Site","Nonstop_Mutation","Translation_Start_Site")) > 0) indel.labels[indel.labels %in% c("Frame_Shift_Del", "Frame_Shift_Ins","Nonsense_Mutation","Splice_Site","Nonstop_Mutation","Translation_Start_Site")] <- "Frame_shift"
                if (sum(indel.labels %in% c("In_Frame_Del", "In_Frame_Ins","Missense_Mutation")) > 0)  indel.labels[indel.labels %in% c("In_Frame_Del", "In_Frame_Ins","Missense_Mutation")] <- "In_frame"
                nonsilent.maf.table[which.indels, type.col] <- indel.labels

                rm(indel.labels)
        }
        rm(type.col)
        rm(class.col)

        Extra.type <- which(!(nonsilent.maf.table[,10] %in% c("SNP", "DNP", "TNP", "ONP", "Frame_shift", "In_frame")))
        if (length(Extra.type) >0){
                for (i in 1:length(Extra.type)){
                        alleles <- c(nonsilent.maf.table[Extra.type[i],11], nonsilent.maf.table[Extra.type[i], 12],nonsilent.maf.table[Extra.type[i], 13])
                        maxchar <- max(nchar(alleles))
                        if ("-" %in% alleles ) {
                                nonsilent.maf.table[Extra.type[i], 10] <- "Frame_shift"
                                if (maxchar %% 3 == 0) nonsilent.maf.table[Extra.type[i], 10] <- "In_frame"
                        } else{
                                nonsilent.maf.table[Extra.type[i], 10] <- "SNP"
                                if (maxchar > 1) nonsilent.maf.table[Extra.type[i], 10] <- "DNP"
                                if (maxchar > 2) nonsilent.maf.table[Extra.type[i], 10] <- "TNP"
                                if (maxchar > 3) nonsilent.maf.table[Extra.type[i], 10] <- "ONP"
                        }
                }
                rm(Extra.type)
                rm(alleles)
                rm(maxchar)
        }

        nonsilent.maf.table = nonsilent.maf.table[!is.na(nonsilent.maf.table[,"Gene"]) & nonsilent.maf.table[,"Chromosome"] %in% c(1:24,"X","Y"),]

        silent.mutation.table <- cbind(Ensembl_gene_id=as.character(silent.maf.table[,"Gene"]), silent.maf.table[,c(5,6,10:13,16)], Protein_position = silent.maf.table[,"Protein_position"], Variant_Classification = silent.maf.table[,"Variant_Classification"])
        nonsilent.mutation.table <- cbind(Ensembl_gene_id=as.character(nonsilent.maf.table[,"Gene"]), nonsilent.maf.table[,c(5,6,10:13,16)], Protein_position = nonsilent.maf.table[,"Protein_position"], Variant_Classification = nonsilent.maf.table[,"Variant_Classification"])



        rm(nonsilent.maf.table)
        rm(silent.maf.table)
        nonsil.maf = maf.dnp.converter(nonsilent.mutation.table)
        sil.maf = maf.dnp.converter(silent.mutation.table)

        rm(nonsilent.mutation.table)
        rm(silent.mutation.table)

        write.table(rbind(cbind(nonsil.maf,1),cbind(sil.maf,2)),file = output,row.names=F,col.names=T,quote=F,sep="\t")
}



maf.dnp.converter <- function(mutab){
        mutab <- as.matrix(mutab)
        mutab.snp= mutab[ mutab[,4]=="SNP" ,]

        if(sum(  mutab[,4]=="In_frame" | mutab[,4]=="Frame_shift"  )>0 )
                mutab.indel=   mutab[ mutab[,4]=="In_frame" | mutab[,4]=="Frame_shift",]   else mutab.indel=NULL

        mutab.dnp= mutab.tnp=mutab.onp=NULL
        if(sum(mutab[,4]=="DNP")>0){
                mutab.dnp=  (as.matrix(mutab[ mutab[,4]=="DNP" ,]))
                if (sum(mutab[,4]=="DNP")==1) {mutab.dnp=  t(as.matrix(mutab[ mutab[,4]=="DNP" ,])) }
                #### convert DNP to SNP ##############
                a= matrix(0, nrow= nrow(mutab.dnp)*2,ncol=ncol(mutab.dnp)   )
                colnames(a)=colnames(mutab)
                a[,1]=  rep(mutab.dnp[,1] ,each=2)
                a[,2]=  rep(mutab.dnp[,2] ,each=2)
                a[,4]=  rep(mutab.dnp[,4] ,each=2)
                a[,8]=  rep(mutab.dnp[,8] ,each=2)
                a[,9]=  rep(mutab.dnp[,9] ,each=2)
                a[,10]=  rep(mutab.dnp[,10] ,each=2)

                a[2*(1:nrow(mutab.dnp))-1,3]=  mutab.dnp[,3]
                a[2*(1:nrow(mutab.dnp)),3]= as.numeric(mutab.dnp[,3]) + 1
                a[,5]= unlist(  strsplit( mutab.dnp[,5]  ,""))
                a[,6]= unlist(  strsplit( mutab.dnp[,6] ,""))
                a[,7]= unlist(  strsplit( mutab.dnp[,7] ,""))
                a = a[(a[,5] == a[,6] & a[,6] == a[,7]) == FALSE,]
                mutab.dnp=a
        }
	if(sum(mutab[,4]=="TNP")>0){
                mutab.tnp=  mutab[ mutab[,4]=="TNP" ,]
                if (sum(mutab[,4]=="TNP")==1) {mutab.tnp=  t(as.matrix(mutab[ mutab[,4]=="TNP" ,])) }
                #### convert TNP to SNP ##############
                a=matrix(0, nrow= nrow(mutab.tnp)*3,ncol=ncol(mutab.tnp)   )
                colnames(a)=colnames(mutab)
                a[,1]=  rep(mutab.tnp[,1] ,each=3)
                a[,2]=  rep(mutab.tnp[,2] ,each=3)
		a[,4]=  rep(mutab.tnp[,4] ,each=3)
                a[,8]=  rep(mutab.tnp[,8] ,each=3)
		a[,9]=  rep(mutab.tnp[,9] ,each=3)
		a[,10]=  rep(mutab.tnp[,10] ,each=3)

                a[3*(1:nrow(mutab.tnp))-2,3]=  mutab.tnp[,3]
                a[3*(1:nrow(mutab.tnp))-1,3]= as.numeric(mutab.tnp[,3]) + 1
                a[3*(1:nrow(mutab.tnp)),3]= as.numeric(mutab.tnp[,3]) + 2
                a[,5]= unlist(  strsplit( mutab.tnp[,5]  ,""))
                a[,6]= unlist(  strsplit( mutab.tnp[,6] ,""))
                a[,7]= unlist(  strsplit( mutab.tnp[,7] ,""))
                a = a[(a[,5] == a[,6] & a[,6] == a[,7]) == FALSE,]
                mutab.tnp=a
        }
        if(sum(mutab[,4]=="ONP")>0)
        {
                mutab.onp.all = mutab[ mutab[,4]=="ONP" ,]
                if (sum(mutab[,4]=="ONP")==1) {mutab.onp.all=  t(as.matrix(mutab[ mutab[,4]=="ONP" ,])) }
                 #### convert ONP to SNP ##############
                maxchar <- apply(cbind(nchar(mutab.onp.all[,5]),nchar(mutab.onp.all[,6]),nchar(mutab.onp.all[,7])),1,max)
                umaxchar <- unique(maxchar)
                b = NULL
                for (k in 1:length(umaxchar)){
                        mutab.onp <- mutab.onp.all[maxchar==umaxchar[k],]
                        if (length(mutab.onp) == 8){mutab.onp = t(as.matrix(mutab.onp))}
                        a= matrix(0, nrow= nrow(mutab.onp)*umaxchar[k],ncol=ncol(mutab.onp)   )
                        colnames(a)=colnames(mutab)
                        a[,1]=  rep(mutab.onp[,1] ,each=umaxchar[k])
                        a[,2]=  rep(mutab.onp[,2] ,each=umaxchar[k])
			a[,4]=  rep(mutab.onp[,4] ,each=umaxchar[k])
                        a[,8]=  rep(mutab.onp[,8] ,each=umaxchar[k])
			a[,9]=  rep(mutab.onp[,9] ,each=umaxchar[k])
			a[,10]=  rep(mutab.onp[,10] ,each=umaxchar[k])

                        for (len in 1:umaxchar[k]){
                                a[umaxchar[k]*(1:nrow(mutab.onp))-(umaxchar[k]-len),3]=  as.numeric(mutab.onp[,3]) + (len-1)
                        }
                        a[,5]= unlist(  strsplit( mutab.onp[,5]  ,""))
                        a[,6]= unlist(  strsplit( mutab.onp[,6] ,""))
                        a[,7]= unlist(  strsplit( mutab.onp[,7] ,""))
                        a = a[(a[,5] == a[,6] & a[,6] == a[,7]) == FALSE,]
                        b <- rbind(b,a)
                }
                mutab.onp <- b
        }

        mutab=rbind(mutab.snp,mutab.dnp,mutab.tnp,mutab.onp,mutab.indel)
        ref=mutab[,5]      ### reference nucleotide
        mut=ref
        temp= (mutab[,6]!=ref)
        mut[temp]=mutab[temp,6]

        temp= (mutab[,7]!=ref)
        mut[temp]=mutab[temp,7]   ### mutated nucleotide

        mutab =  cbind(mutab,mut)
        colnames(mutab)[11] = "Mut"

        return(mutab)
}

maf.cleanup(maf.file,output)
