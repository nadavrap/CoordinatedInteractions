rm( list=ls() )
load('Rdata/setup.Rdata')
source( 'simfxn.R' )

iter		<- as.numeric( commandArgs(TRUE)[[1]] )
mode		<- as.numeric( commandArgs(TRUE)[[2]] )

savefile	<- paste0( 'Rdata/', mode, '_', iter, '.Rdata')
sinkfile	<- paste0( 'Rout/' , mode, '_', iter, '.Rout' )
if( file.exists(sinkfile) | file.exists(savefile) ) stop()
sink( sinkfile )

out	<- rbind(
	main(seed=iter, N=Ns[mode], M, simtype=simtypes[mode], h2=h2s[mode], gamma=gams[mode], maftype=maftypes[mode], adjust_pcs=adj_pcs[mode], oracle=oracle[mode], trips=trips[mode], epitype='add'),
	main(seed=iter, N=Ns[mode], M, simtype=simtypes[mode], h2=h2s[mode], gamma=gams[mode], maftype=maftypes[mode], adjust_pcs=adj_pcs[mode], oracle=oracle[mode], trips=trips[mode], epitype='isotropic'),
	main(seed=iter, N=Ns[mode], M, simtype=simtypes[mode], h2=h2s[mode], gamma=gams[mode], maftype=maftypes[mode], adjust_pcs=adj_pcs[mode], oracle=oracle[mode], trips=trips[mode], epitype='ci')
)
print( out )

save( out, file=savefile )
sink()
