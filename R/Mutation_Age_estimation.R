

#Modified from Gandolfo, L. C., Bahlo, M., and Speed, T. P. (2014).
#Dating rare mutations from small samples with dense marker data. Genetics, 197(4): 1315-27. 


mutation_age_estimation = function(segments = segments,
                                   confidence.coefficient = 0.95,
                                   chance.sharing.correction = T,
                                   median.allele.frequency = 0.5,
                                   markers.on.chromosome = 23144,
                                   length.of.chromosome = 278.1,
                                   focus = 156.107442,
                                   patient_name = 'patient1',
                                   map_file = 'genetic_map_chr1_b36.txt'){

  
#physical to genetic distance
segments = physical_genetic(segments,map_file=map_file,focus=focus)


#Chance sharing correction
if (chance.sharing.correction == TRUE){
	e = 0.01
	p = (median.allele.frequency)^2 + (1-median.allele.frequency)^2
	phi = (length.of.chromosome/100)/markers.on.chromosome
	loci = log(e)/log(p)
	cs.correction = loci*phi
	}
if (chance.sharing.correction == FALSE){cs.correction = 0}

#Age estimation and confidence intervals
cc = confidence.coefficient
l.lengths = (1/100)*segments$l.lengths 
r.lengths = (1/100)*segments$r.lengths
n = length(l.lengths)

#Assuming an 'independent' genealogy
if (n < 10){i.cs.correction = 0}
if (n >= 10){i.cs.correction = cs.correction}

length.correction = (sum(l.lengths) + sum(r.lengths) - 2*(n-1)*i.cs.correction)/(2*n)

sum.lengths = sum(l.lengths) + sum(r.lengths) + 2*length.correction - 2*(n-1)*i.cs.correction
b.c = (2*n-1)/(2*n)
i.tau.hat <- (b.c*2*n)/sum.lengths
g_l <- qgamma(shape=2*n,scale=1/(2*n*b.c),((1-cc)/2))
g_u <- qgamma(shape=2*n,scale=1/(2*n*b.c),(cc+(1-cc)/2))		
i.l = g_l*i.tau.hat
i.u = g_u*i.tau.hat

#Assuming a 'correlated' genealogy
length.correction = (sum(l.lengths) + sum(r.lengths) - 2*(n-1)*cs.correction)/(2*n)

longest.l.lengths = match(sort(l.lengths,decreasing=TRUE)[1], l.lengths)
l.lengths[longest.l.lengths] <- l.lengths[longest.l.lengths] + length.correction + cs.correction
longest.r.lengths = match(sort(r.lengths,decreasing=TRUE)[1], r.lengths)
r.lengths[longest.r.lengths] <- r.lengths[longest.r.lengths] + length.correction + cs.correction

lengths = l.lengths + r.lengths
lengths = lengths - 2*cs.correction
rho.hat = (n*(mean(lengths))^2 - var(lengths)*(1+2*n))/(n*(mean(lengths))^2 + var(lengths)*(n-1))
n.star = n/(1+(n-1)*rho.hat)
	if (n.star > n) {n.star = n}
	if (n.star < -n) {n.star = -n}		
b.c = (2*n.star-1)/(2*n.star)
c.tau.hat = (b.c*2*n)/sum(lengths)
	if (rho.hat < -2/(n-1)){n.star = n/(1+(n-1)*abs(rho.hat))}
	if (-2/(n-1) <= rho.hat & rho.hat < -1/(n-1)){n.star = n}
g_l = qgamma(shape=2*n.star,scale=1/(2*n.star*b.c),(1-cc)/2)
g_u = qgamma(shape=2*n.star,scale=1/(2*n.star*b.c),cc+(1-cc)/2)
c.l = g_l*c.tau.hat
c.u = g_u*c.tau.hat

#Print results
#print(paste("Assuming an 'independent' genealogy: age estimate =", round(i.tau.hat, digits=1) ,"generations, CI: ", paste("(",round(i.l, digits=1),",",round(i.u, digits=1),")",sep="")), quote=FALSE) 
#print(paste("Assuming a 'correlated' genealogy: age estimate =", round(c.tau.hat, digits=1) ,"generations, CI: ", paste("(",round(c.l, digits=1),",",round(c.u, digits=1),")",sep="")), quote=FALSE)

out = data.frame(patient = c(patient_name,patient_name), genealogy = c('independent','correlated'), generations = c(i.tau.hat,c.tau.hat), CI_minus = c(i.l,c.l), CI_plus = c(i.u,c.u))
                 
return(out)
}



#Calculate genetic distance from physical distance based on a genetic map from here: https://ftp.ncbi.nlm.nih.gov/hapmap/recombination/latest/rates/genetic_map_chr1_b36.txt
physical_genetic = function(segments = 'this',focus = 156107442,  map_file = 'genetic_map_chr1_b36.txt'){
  
  map = read.table(map_file,header = T)
  map$position = map$position/1000000
    
  segments$l.lengths = 0
  segments$r.lengths = 0
  
  focal_cM = map[map$position > focus ,]
  focal_cM = focal_cM[1,3]
  
  for(i in 1:nrow(segments)){
    start = map[map$position > segments$V1[i] ,]
    start = ifelse(nrow(start)==0,map[1,3],start[1,3])
    
    end = map[map$position > segments$V2[i] ,]
    end = ifelse(nrow(end)==0,map[nrow(map),3],end[1,3])
    
    segments$l.lengths[i] = focal_cM - start
    segments$r.lengths[i] = end - focal_cM 
  }
  
  return(segments)
}

