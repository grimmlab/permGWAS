### R function to plot GWAS results from  permGWAS ##

# plot_permGWAS(output=read.csv('p_values_phenotype_value.csv'),q2=read.csv('min_p_values_phenotype_value.csv'))



plot_permGWAS<-function(output,q2,h=3,black=TRUE,colSet='Dark2',mafi=0.01,thres=NA,max.y=NA,lower.limit=0.01,title=TRUE,name=NA) {

  colnames(output)[h]<-'Plot'
  
  q2<- q2[order(q2[,2]),]  
  q2<- -log10 (q2[ceiling(0.05*nrow(q2)),2])
  
  output_<-subset(output,output$maf>mafi)
  new_<-subset(output_,output_$Plot<lower.limit)
  if (is.na(thres)==TRUE) {thres<--log10(0.05/nrow(output_))}
  
  
  
  
  if(is.na(max.y)==TRUE) {
    max.y<-ceiling(-log10(min(new_[,h]))+1)
    if(max.y< q2)
    {max.y=q2}} 
  
  output_<-new_[order(new_$POS),]
  output_ok<-output_[order(output_$CHR),]
  
  maxpos<-c(0,cumsum(aggregate(output_ok$POS,list(output_ok$CHR),max)$x+max(cumsum(aggregate(output_ok$POS,list(output_ok$CHR),max)$x))/100))
  if(black==TRUE) {
    plot_col<-rep(c('gray10','gray60'),ceiling(max(unique(output_ok$CHR))/2))
  } else { 
    require(RColorBrewer)
    colorCount = length(unique(output_ok$CHR))
    getPalette = colorRampPalette(brewer.pal(colorCount,colSet))
    
    plot_col<-getPalette(colorCount)}
  
  size<-aggregate(output_ok$POS,list(output_ok$CHR),length)$x
  
  a<-rep(maxpos[1],size[1])
  b<-rep(plot_col[1],size[1])
  for (i in 2:max(unique(output_ok$CHR))){
    a<-c(a,rep(maxpos[i],size[i]))
    b<-c(b,rep(plot_col[i],size[i]))}
  
  output_ok$xpos<-output_ok$POS+a
  output_ok$col<-b
  
  d<-(aggregate(output_ok$xpos,list(output_ok$CHR),min)$x+aggregate(output_ok$xpos,list(output_ok$CHR),max)$x)/2
  
  plot(output_ok$xpos,-log10(output_ok$Plot),col=output_ok$col,pch=16,ylab='-log10(pval)',xaxt='n',xlab='chromosome',axes=FALSE,cex=1.2,ylim=c(-log10(lower.limit),max.y))
  axis(1,tick=FALSE,at=d,labels=c(1:max(unique(output_ok$CHR))))
  axis(2,lwd=2)
  abline(h=thres,lty=9,col='#e41a1c',lwd=2)
  abline(h=q2,lty=2,col='#377eb8',lwd=2)
  
  if(title==TRUE) {
    title(main=name)}

}

