#this package assumes a linux operating system with the following operational programs:
#mutt
#cvlc (ie. a terminal version of vlc)
#att cell phone servace
stamp = proc.time()

ding = function(num){
	system(sprintf('cvlc DING%d.mp3 &', num))
}

txtMsg = function(output, outName, number){
	#output should be a string version of all of the information that you want to write to file with \n characters separating lines and \escaped special characters
        #outName should be the path to an output filename as a string
        #number should be given as a string
	tempo = proc.time()-stamp	

	sink(outName)
	writeLines(output)
	writeLines('\n')
	print(tempo)
	sink()	

	emailAddress = sprintf("%s@txt.att.net", number)
	system(sprintf('echo "" | mutt -s "R Update: %s" -a %s %s', date(), outName, emailAddress))
}

email = function(output, outName, emailAddress){
	#output should be a string version of all of the information that you want to write to file with \n characters separating lines and \escaped special characters
	#outName should be the path to an output filename as a string
	#emailAddress should be given as a string
	tempo = proc.time()-stamp

	sink(outName)
	writeLines(output)
	writeLines('\n')
	print(tempo)
	sink()

	system(sprintf('echo "" | mutt -s "R Update: %s" -a %s %s', date(), outName, emailAddress))
}
