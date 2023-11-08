#!/usr/bin/ruby

class Vec
	def initialize(name, size)
		@name = name
		@size = size
		@i = Array.new(size)
		@v = Array.new(size)
		@k = 0
	end

	def add(i, v)
		@i[@k] = i
		@v[@k] = v
		@k += 1
	end

	def printVec
		print "#{@name} <- array(0, #{@size})\n"
		print "i <- c("
		print @i.join(",")
		print ")\n"
		print "x <- c("
		print @v.join(",")
		print ")\n"
		print "#{@name}[i] <- x\n"
	end
end

res = []
while line = gets do
	if line =~ /^#/
		a = line.split(/\s+/)
		name = a[1]
		line = gets
		a = line.split(/\s+/)
		size = a[2].to_i
		m = Vec.new(name, size)
		res << m
	else
		a = line.split(/\s+/)
		m.add(a[1], a[2])
	end
end

res.each do |x|
	x.printVec
end

