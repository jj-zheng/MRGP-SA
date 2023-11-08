#!/usr/bin/ruby

class Mat
	def initialize(name, x, y, nnz)
		@name = name
		@x = x
		@y = y
		@nnz = nnz
		@i = Array.new(nnz)
		@j = Array.new(nnz)
		@v = Array.new(nnz)
		@k = 0
	end

	def add(i, j, v)
		@i[@k] = i
		@j[@k] = j
		@v[@k] = v
		@k += 1
	end

	def printMat
		print "i <- c("
		print @i.join(",")
		print ")\n"
		print "j <- c("
		print @j.join(",")
		print ")\n"
		print "x <- c("
		print @v.join(",")
		print ")\n"
		print "#{@name} <- sparseMatrix(i=i, j=j, x=x, dims=c(#{@x}, #{@y}))\n"
	end
end

res = []
while line = gets do
	if line =~ /^#/
		a = line.split(/\s+/)
		name = a[1]
		line = gets
		a = line.split(/\s+/)
		x = a[2].to_i
		y = a[3].to_i
		nnz = a[4].to_i
		m = Mat.new(name, x, y, nnz)
		res << m
	else
		a = line.split(/\s+/)
		m.add(a[1], a[2], a[3])
	end
end

res.each do |x|
	x.printMat
end

