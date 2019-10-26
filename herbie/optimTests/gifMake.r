gifMake = function(itMeat, itRange, name, delay, width=700, height=700, max=999){
	#itMeat: a function that you want to iterate over the it range, your plotting functions should occure here.
	#itRange: a vector to iterate itMeat over.
	#name: an output name given as a string
	#delay: the delay between frames given in seconds.
	
	count=(max+1)*9
	if (length(itRange)>max) {
		warning('Number of requested iterations exceeds maximum.\nDecrease length(itRange), or Increase max.')
		break
	}
	for (it in itRange){
		jpeg(sprintf('./%s/%4s.jpeg', name, count), width=width, height=height, quality=100)
		itMeat()
		dev.off()
		count = count+1
	}
	
	#make the gif in the current directory for convenience
	system(sprintf('convert -delay %d00 ./%s/*.jpeg ./%s.gif', delay, name, name))
	#also put the gif in the storage directory for goood form
	system(sprintf('cp ./%s.gif ./%s/%s.gif', name, name, name) )
}
