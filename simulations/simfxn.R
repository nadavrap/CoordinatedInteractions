getmafs	<- function(maftype,M,tmafs,Fst){
	if(maftype=='unif'){
		rawmafs	<- rbeta(M, 1, 1)
	} else if(maftype == 'strat'){
		rawmafs	<- sapply( tmafs, function(tmaf) rbeta(1,tmaf*(1-Fst)/Fst,(1-tmaf)*(1-Fst)/Fst) )
	} else if(maftype == 'UKBB'){
		rawmafs <- sample( read.table("./ukbb_imp_MAF0.001.afreq.gz", header = F, stringsAsFactors=F)[,5], M )
	}
	sapply( rawmafs, function(x) max( min(x,.99), .01 ) )
}

genofxn	<- function(simtype='simple',maftype='unif',Fst,N,M,...){
	tmafs <- getmafs(maftype,M)
	if( simtype=='simple' ){
		genos	<- sapply( tmafs, function(x) rbinom(N,2,x) )
	} else if( simtype=='am' ){
		genos	<- sapply( tmafs, function(x) rbinom(2*N,2,x) ) # parent genotypes
		pheno	<- phenofxn( 2*N, M, scale(genos), parents=TRUE, ... )
		### ------ Mate parents by phenotype ------ ###
		sort.genos <- genos[order(pheno), ] # sort genos by sorted phenotypes
		dgenos <- sort.genos[seq(1,2*N-1,by=2), ] # sires
		mgenos <- sort.genos[seq(2,2*N  ,by=2), ] # dams
		### ------- Make child genotypes based on parent genotypes ------ ###
		mallele <- matrix( rbinom(N*M,1,mgenos/2), N, M )
		dallele <- matrix( rbinom(N*M,1,dgenos/2), N, M )
		genos		<- mallele + dallele # child genotypes
	} else if( simtype=='strat' ){
		tmafs1 <- getmafs('strat',M,tmafs,Fst)
		tmafs2 <- getmafs('strat',M,tmafs,Fst)
		genos1<- t(sapply(1:(N/2),function(x) rbinom(M,2,tmafs1) ))
		genos2<- t(sapply(1:(N/2),function(x) rbinom(M,2,tmafs2) ))
		genos	<- rbind(genos1,genos2)
	}
	scale(genos)
}

phenofxn	<- function( seed, N, M, genos, epitype, h2, h2x=.1, gamma, betaz, trips=FALSE, propInt=.01, Mpath=M*3/10, propInt3=0.0001, z, parents=FALSE ){
	set.seed(seed)
	betag <- sqrt(h2/M) * rnorm(M)
	addg	<- genos%*%betag
	if( epitype == 'add' ){
		intrxn	<- 0
		sig2e		<- 1-h2-betaz^2
	} else if( epitype == 'isotropic' & !trips ){
		Mx	<- M^2*propInt
		Gdoubs	<- sapply( 1:Mx, function(j){
			doub	<- sample( M, 2, replace=F )
			genos[,doub[1]] * genos[,doub[2]]
		})
		betasx	<- sqrt(h2x/Mx) * rnorm(Mx)
		intrxn	<- Gdoubs%*%betasx
		sig2e		<- 1-h2-h2x-betaz^2
	} else if( epitype == 'isotropic' & trips ){
		Mx	<- M^3*propInt3
		Gtrips	<- sapply( 1:Mx, function(j){
			trip	<- sample( M, 3, replace=F )
			genos[,trip[1]] * genos[,trip[2]] * genos[,trip[3]]
		})
		betasx	<- sqrt(h2x/Mx)*rnorm(Mx)
		intrxn	<- Gtrips %*% betasx
		sig2e		<- 1-h2-h2x-betaz^2
	} else if( epitype == 'ci' & !trips ){
		doubs	<- 1:(2*Mpath)
		paths	<- sapply( 1:2, function(k){
			doubloc	<- 1:Mpath +(k-1)*Mpath
			genos[,doubloc]%*%betag[doubloc]
		})
		intrxn<- gamma * scale(apply( paths, 1, prod ))
		sig2e	<- 1-h2-gamma^2-betaz^2
	} else if( epitype == 'ci' & trips ){
		trips	<- sample( M, 3*Mpath, rep=F )
		paths	<- sapply( 1:3, function(k){
			triploc	<- 1:Mpath +(k-1)*Mpath
			genos[,triploc]%*%betag[triploc]
		})
		intrxn<- gamma * scale(apply( paths, 1, prod ))
		sig2e	<- 1-h2-gamma^2-betaz^2
	}
	#print( c( mean( addg^2 ), mean( intrxn^2 ), mean( (z*betaz)^2 ), mean( (sqrt(sig2e)*rnorm(N) )^2 ) ))
	if( parents ){
		set.seed( seed + 1e6 )	### to ensure different epsilons,
		z	<- rep(z,2)						### because #parents = 2N
	}
	as.numeric( addg + intrxn + z*betaz + sqrt(sig2e)*rnorm(N) )
}

testfxn	<- function(pheno,genos,adjust_pcs,oracle){
	#select half of the individuals to run a marginal gwas and the other half to test for interaction
	N			<- length(pheno)
	test	<- sample(1:N, N/2)
	train <- setdiff(1:N,test)

	# run a gwas on the train data
	if(!adjust_pcs){
		betas	<- sapply( 1:M, function(k) summary(lm(pheno[train]~genos[train,k]))$coef[2,1] )
	} else {
		pc		<- svd(genos[train,])$u[,1]
		betas	<- sapply( 1:M, function(k) summary(lm(pheno[train]~genos[train,k]+pc))$coef[2,1] )
	}

	# compute PRS and even/odd versions
	genos	<- genos			[test,]
	pheno	<- scale(pheno[test])
	if( oracle ){
		Mpath	<- 300
		even	<- 1:Mpath
		odd		<- Mpath+1:Mpath
	} else {
		even	<- sample( M, M/2 )
		odd		<- setdiff( 1:M, even )
	}
	prs		<- scale( genos       %*%betas			)
	prse	<- scale( genos[,even]%*%betas[even])
	prso	<- scale( genos[,odd ]%*%betas[odd ])
	eo		<- prse*prso 

	# test out of sample: PRS effect; EO corr; EO Test
	if(!adjust_pcs){
		m.prs	<- summary(	lm( pheno	~ prs))
		m.cor <- summary(	lm( prso	~ prse))
		m.int	<- summary(	lm( pheno	~ 1+eo+prse+prso))
	} else {
		pc		<- svd(genos)$u[,1]
		m.prs	<- summary(	lm( pheno	~ prs+pc))
		m.cor <- summary(	lm( prso	~ prse+pc))
		m.int	<- summary(	lm( pheno	~ 1+eo+prse+prso+pc))
	}

	return(c(
		m.prs$coefficients[2,3:4],#PRS t, p-value
		m.cor$coefficients[2,3:4],#eo-corr t, p-value
		m.int$coefficients[2,3:4] #eo-interxn t, p-value
	))
}

main <- function( seed, N, M, simtype='simple', maftype='unif', Fst=0.1, adjust_pcs=FALSE, betaz=ifelse( simtype=='strat', .2, 0 ), oracle, ... ){
	set.seed(seed)
	z			<- scale( rep(0:1,each=N/2) )
	genos	<- genofxn(simtype,maftype,Fst,N,M,z=z,betaz=betaz,seed=seed,...)
	pheno	<- phenofxn( seed, N, M, genos, z=z,betaz=betaz,... )
	testfxn( pheno, genos, adjust_pcs, oracle )
}
