rm(list=ls())

data=read.csv(file="grades.csv",head=TRUE,sep=",")

dev.new()
plot(density(data$Midterm..100.),
	xlab='Midterm Score',
	main='AMS 7 Spring 2013\n Midterm Results')

dev.new()
hist(data$Midterm..100., breaks=20,
	xlab='Midterm Score',
	main='AMS 7 Spring 2013\n Midterm Results')