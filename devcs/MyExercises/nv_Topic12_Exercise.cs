using System;


namespace nvolfango_Topic12_Exercises
{
	// Monte Carlo Integration
	class MonteCarloIntegration
	{
		public static void MainProgram()
		{
			// Initialise inputs
			int n = 100000;		// Number of samples to be evaluated
			int upper = 10;		// Upper limit of integral
			int lower = 0;		// Lower limit of integral
			double mc_error;	// Estimate of error

			double integral;
			double true_value;

			// MC Integration for first function:
			GeneralFunction f1 = new GeneralFunction((double x) => Math.Sqrt(x + Math.Sqrt(x)));
			integral = EvaluateMonteCarloIntegral(f1, out mc_error, lower, upper, n);
			true_value = 25.52664964560352;

			Console.WriteLine("Function 1:");
			Console.WriteLine("Monte Carlo Integral = {0}", integral);
			Console.WriteLine("Error = {0}", mc_error);
			Console.WriteLine("True value of integral = {0}\n", true_value);

			// MC Integration for second function:
			GeneralFunction f2 = new GeneralFunction((double x) => Math.Cos(x + Math.Sin(x)));
			integral = EvaluateMonteCarloIntegral(f2, out mc_error, lower, upper, n);
			true_value = -4.722560808683198;

			Console.WriteLine("Function 2:");
			Console.WriteLine("Monte Carlo Integral = {0}", integral);
			Console.WriteLine("Error = {0}", mc_error);
			Console.WriteLine("True value of integral = {0}", true_value);
		}


		delegate double GeneralFunction(double x);


		static double EvaluateMonteCarloIntegral(GeneralFunction function, out double error, int lower_limit, int upper_limit, int n)
		// Calculates the Monte Carlo integral for a function (input and output must be a double)
		// upper_limit: upper bound of the integral
		// lower_limit: lower bound of the integral
		// n: Number of (random) points to be evaluated
		
		{
			double sum = 0;
			double sum_sq = 0;		// Used for calculatnig error
			double f_val;
			double f_bar;
			double f_bar_sq;		// Used for calculating error
			double integral;

			Random r = new Random();

			for (int i = 0; i < n; i++)
			{
				f_val = (r.NextDouble() * (upper_limit - lower_limit)) + lower_limit;
				sum += function(f_val);
				sum_sq += function(f_val) * function(f_val);
			}
			f_bar = sum / n;
			f_bar_sq = sum_sq / n;

			integral = (upper_limit - lower_limit) * f_bar;
			error = (upper_limit - lower_limit) * Math.Sqrt((f_bar_sq - f_bar * f_bar) / n);

			return integral;
		}
	}
}