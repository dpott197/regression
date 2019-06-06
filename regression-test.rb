require 'matrix'
@factors = Matrix[[1,40,10]]
@cofactors = Matrix[[5.43405,-0.085930,-0.118865],[-0.085930,0.00147070,0.00165094],[-0.118856,0.00165094,0.00359276]]

load "./regression.rb"

y  = [12.4, 11.7, 12.4, 10.8,  9.4,  9.5,  8.0,  7.5]
x1 = [28.0, 28.0, 32.5, 39.0, 45.9, 57.8, 58.1, 62.5]
x2 = [  18,   14,   24,   22,    8,   16,    1,    0]

xs = y.count.times.collect{|i| [1,x1[i],x2[i]]}
ys = y.each.collect{|y| [y]}

puts @regression = Regression.new(xs,ys)

puts @regression.degrees_of_freedom
# => 5

puts @regression.estimate(Matrix[[1,40,10]])
# => 10.333132089815114

puts @regression.distance_value(Matrix[[1,40,10]])
# => 0.2156687027981243

puts @regression.normalize(Matrix[[1,40,10]],10.333132089815114)
# =>  0.5

puts @regression.normalize(Matrix[[1,40,10]],9.293)
# => 0.975

puts @regression.normalize(Matrix[[1,40,10]],11.374)
# => 0.025

puts @regression.limit(Matrix[[1,40,10]])
#=> 9.29276762699449

puts @regression.limit(Matrix[[1,40,10]],0.975)
#=> 11.37349655263574
