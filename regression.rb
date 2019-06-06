class Regression

  require 'distribution'
  require 'matrix'

  attr_accessor *[
    :adjusted_multiple_coefficient_of_determination,
    :coefficients,
    :covariance_matrix,
    :degrees_of_freedom,
    :f_statistic,
    :errors,
    :explained_variation,
    :multiple_coefficient_of_determination,
    :standard_error,
    :sum_of_squared_errors,
    :t_values,
    :total_variation,
    :unexplained_variation,
    :x_averages,
    :x_matrix,
    :xs,
    :y_average,
    :y_estimates,
    :y_matrix,
    :ys,
  ]

  def initialize(xs,ys)
    if xs.count == ys.count
      run(xs,ys)
    else
      puts "Error: The row counts of x and y do not match."
    end
  end

  def mismatch
    mismatch = false
    @xs.each do |x1|
      @xs.each do |x2|
        if x1.size != x2.size
          mismatch = true
        end
      end
    end
    mismatch
  end

  def run(xs,ys)
    to_dataset(xs,ys)
    if data_are_acceptable?
      coefficients
       y_estimates
      errors
      sum_of_squared_errors
      y_average
      total_variation
      explained_variation
      unexplained_variation
    end
  end

  def coefficients
    @coefficients ||= covariance_matrix.inverse*@x_matrix.transpose*@y_matrix
  end

  def x_matrix(xs=nil)
    @x_matrix ||= Matrix[*xs]
  end

  def y_matrix(ys=nil)
    @y_matrix ||= Matrix[*ys]
  end

  def covariance_matrix
    @covariance_matrix ||= @x_matrix.transpose*@x_matrix
  end

  def data_are_acceptable?
    covariance_matrix.regular? && degrees_of_freedom.to_f > 0.0
  end

  def to_dataset(xs,ys)
    @xs = xs; @ys = ys.flatten; x_matrix(xs); y_matrix(ys);
  end

  def degrees_of_freedom
    @degrees_of_freedom ||= (@x_matrix.row_size-@x_matrix.column_size).to_f
  end

  def estimate(xs)
    (Matrix[*xs]*coefficients)[0,0].to_f
  end

  def forecast(xs)
    estimate(xs)
  end

  def predict(xs)
    estimate(xs)
  end

  def y_estimates
    @y_estimates ||= @xs.collect{|x| estimate([x])}
  end

  def estimates
    y_estimates
  end

  def errors
    @errors ||= @y_estimates.count.times.collect{|i| @ys[i].to_f-@y_estimates[i].to_f}
  end

  def sum_of_squared_errors
    @sum_of_squared_errors ||= @errors.inject(0){|sum,error| sum.to_f+error.to_f**2.0}
  end

  def y_average
    @y_average ||= @ys.flatten.inject(0){|sum,actual| sum+actual.to_f}/@ys.flatten.count
  end

  def total_variation
    @total_variation ||= @ys.inject(0){|sum,y| sum + (y.to_f-y_average.to_f)**(2.0)}
  end

  def explained_variation
    @explained_variation ||= @y_estimates.inject(0){|sum,y_estimate| sum + (y_estimate.to_f-y_average.to_f)**(2.0)}
  end

  def unexplained_variation
    @unexplained_variation ||= sum_of_squared_errors
  end

  def standard_error
    @standard_error ||= Math.sqrt(unexplained_variation.to_f/degrees_of_freedom.to_f)
  end

  def multiple_coefficient_of_determination
    @multiple_coefficient_of_determination ||= explained_variation.to_f/total_variation.to_f
  end

  def multiple_correlation_coefficient
    @multiple_correlation_coefficient ||= Math.sqrt(multiple_coefficient_of_determination)
  end

  def r_squared
    multiple_coefficient_of_determination
  end

  def r2
    multiple_coefficient_of_determination
  end

  def adj_r2
    adjusted_multiple_coefficient_of_determination
  end

  def r
    multiple_correlation_coefficient
  end

  def adjusted_r_squared
    adjusted_multiple_coefficient_of_determination
  end

  def f_statistic
    numerator = explained_variation/k
    denominator = unexplained_variation/(n - (k+1.0))
    numerator/denominator
  end

  def adjusted_multiple_coefficient_of_determination
    r2 - (1.0-r2.to_f)*((k.to_f)/(n.to_f-k.to_f-1.0))
  end

  def n
    @x_matrix.row_size.to_f
  end

  def k
    @x_matrix.column_size.to_f-1
  end

  def df
    degrees_of_freedom
  end

  def f
    f_statistic
  end

  def f_stat
    f_statistic
  end

  def distance_value(xs)
    (Matrix[*xs]*covariance_matrix.inverse*Matrix[*xs].transpose)[0,0]
  end

  def penalty_factor(xs,prediction_penalty=1)
    Math.sqrt(prediction_penalty+distance_value(xs))
  end

  def penalized_standard_error(xs,prediction_penalty)
    standard_error*penalty_factor(xs,prediction_penalty)
  end

  def half_width(xs,probability,prediction_penalty)
    Distribution::T.p_value(probability,degrees_of_freedom)*penalized_standard_error(xs,prediction_penalty)
  end

  def limit(xs,probability=0.025,prediction_penalty=1)
    if coefficients.nil?
      run
    end
    (Matrix[*xs]*coefficients)[0,0].to_f+half_width(xs,probability,prediction_penalty)
  end

  def normalize(xs,y,prediction_penalty=1)
    if @coefficients.nil?
      run
    end
    Distribution::T.cdf(((Matrix[*xs]*@coefficients)[0,0].to_f-y)/(penalized_standard_error(xs,prediction_penalty)),degrees_of_freedom)
  end

  def probability(xs,y,prediction_penalty=1)
    normalize(xs,y,prediction_penalty)
  end

  def t_values
    @t_values ||= covariance_matrix.row_size.times.collect{|i| coefficients[i,0].to_f/(standard_error.to_f*Math.sqrt(covariance_matrix.inverse[i,i].to_f))}
  end

  def t_statistics
    t_values
  end

  def t_stats
    t_values
  end

  def significant?(percentage=0.95,column=nil)
    if column.nil?
      Distribution::F.cdf(f,k,df) > percentage
    else
      Distribution::T.cdf(t_values[column],degrees_of_freedom) > percentage
    end
  end

end
