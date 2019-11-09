using System;

namespace nv_tutorial03
{
	class nv_tutorial03_exercise01
	{
		public static void MainProgram()
		{
			Random rand = new Random();
			const int n = 10000;
			double[] my_array = new double[n];
			double upper = 50;
			double lower = -50;

			for (int i = 0; i < my_array.Length; i++)
			{
				my_array[i] = rand.NextDouble() * (upper - lower) + lower;
			}

			double mean = Mean(my_array);
			double sd = SD(my_array);
			Console.WriteLine("Mean: {0}, Standard Deviation: {1}", mean, sd);
		}

		public static double Mean(double[] array)
		{
			double sum = 0;
			double mean;

			for (int i = 0; i < array.Length; i++)
			{
				sum += array[i];
			}

			mean = sum / array.Length;

			return mean;
		}

		public static double SD(double[] array)
		{
			double mean = Mean(array);
			double diff;
			double sq_diff;
			double sum = 0;
			double sd;

			for (int i = 0; i < array.Length; i++)
			{
				diff = array[i] - mean;
				sq_diff = diff * diff;
				sum += sq_diff;
			}

			sd = sum / array.Length;

			return Math.Sqrt(sd);
		}
	}

	class nv_tutorial03_exercise02
	{
		// Already done in tutorial05 (Matrix class)
	}
}
