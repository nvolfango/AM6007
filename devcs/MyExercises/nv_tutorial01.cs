/* Each exercise is enclosed within a class, and the entry function called "MainProgram()" for each exercise,
 * and can be executed with the namespace nv_tutorial01, i.e. "using nv_tutorial01", as follows:
 * public static void Main(string[] args)
 * {
 *		nv_tutorial01_exercise01.MainProgram();
 * }
 */

using System;
using System.Collections.Generic;
using System.Text;

namespace nv_tutorial01
{
	class nv_tutorial01_exercise01
	// Calculator
	{
		public static void MainProgram()
		{
			float a, b = 0;
			string num_input;   // Stores user number input

			// Ask user for the two numbers to input
			Console.WriteLine("Please give me two numbers: ");

			Console.Write("First number: ");
			num_input = Console.ReadLine();
			a = float.Parse(num_input); ;

			Console.Write("Second number: ");
			num_input = Console.ReadLine();
			b = float.Parse(num_input);

			// Ask user for the operation they would like to use with the two numbers.
			Console.WriteLine("What operation would you like to run on the two numbers?");
			Console.WriteLine("Choices: ");

			// List out all possible operations.
			List<string> operationList = new List<string>() { "addition", "subtraction", "division", "multiplication" };
			operationList.ForEach(op => Console.WriteLine("\t{0}", op));

			bool valid_input = false;   // Used to ensure that user enters a valid operation.
			float answer = 0;           // Final answer from chosen operation

			// Executes the user's chosen operation if input is valid. Otherwise, keep requesting another input.
			while (!valid_input)
			{
				string operation = Console.ReadLine();
				switch (operation)
				{
					case "addition":
						answer = a + b;
						valid_input = true;
						break;
					case "subtraction":
						answer = a - b;
						valid_input = true;
						break;
					case "division":
						answer = a / b;
						valid_input = true;
						break;
					case "multiplication":
						answer = a * b;
						valid_input = true;
						break;
					default:
						Console.WriteLine("Invalid input. Please choose an operation in the list of choices.");
						break;
				}
			}
			Console.WriteLine("Answer: {0}", answer);
		}
	}

	class nv_tutorial01_exercise02
	// Solve transcendental equations
	{
		public static void MainProgram()
		{
			const int max_iterations = 1000;    // Maximum allowable number of iterations
			const double eps = 1E-3;            // Tolerance for acceptance of numerical solution
			double answer = 0;                  // Final solution
			int num_iters;

			// Initial parameter
			double x = 0.53;                    // Initial guess for x.
			double y = 0.5;

			for (num_iters = 0; num_iters <= max_iterations; num_iters++)
			{
				y = Equation(x);           // Function evaluated at initial guess.
				double d = Distance(x, y);
				if (d <= eps)
				{
					answer = x;
					break;
				}
				else
				{
					x = (x + y) / 2;
				}
			}

			if (num_iters == max_iterations + 1)
			{
				Console.WriteLine("Maximum iterations reached. Could not find a solution.");
			}
			else
			{
				Console.WriteLine("Final answer: {0:N3}", answer);
				Console.WriteLine("Number of iterations: {0}", num_iters);
			}
		}

		static double Equation(double x)
		// Equation for which we want to solve x.
		{
			//return Math.Cos(x);
			//return (Math.Exp(x)-2);
			return (2 * Math.Cos(x));
		}

		static double Distance(double y1, double y2)
		{
			return Math.Abs(y1 - y2);
		}
	}

	class nv_tutorial01_exercise03
	// Infinite Monkey Theorem
	{
		public static void MainProgram()
		{
			// Ask user for the desired length of string.
			Console.Write("Desired string length: ");
			int string_length = int.Parse(Console.ReadLine());

			// Ask user for name they want to find in the string.
			string name;
			Console.Write("Name to search for within string: ");
			name = Console.ReadLine();

			// Generate the random character sequence.
			string monkey_string = GenerateRandomString(string_length);

			// Check whether or not user name is in the sequence.
			if (monkey_string.IndexOf(name) == -1)
			{
				Console.WriteLine("Name \"{0}\" not found in monkey string.\n", name);
			}
			else
			{
				int name_index = monkey_string.IndexOf(name);
				Console.WriteLine("Your name (\"{0}\") was first found in the position no. {1}!!!", name, name_index + 1);
				Console.WriteLine("\"..." + monkey_string.Substring(name_index - 3, name.Length + 6) + "...\"\n");
			}

			// Calculate probability of finding name in the random string.
			double p = CalculateMonkeyProbability(string_length, name.Length);
			Console.WriteLine("Probability of finding \"{0}\" in a sequence of {1} characters: {2}", name, string_length, p);
			Console.WriteLine("That is, you have a 1 in {0} chance of finding \"{1}\" in the string.", Math.Ceiling(1 / p), name);
		}

		static string GenerateRandomString(int string_length)
		{
			const byte start = 97;              // lowercase letters
			const byte end = 123;               // start from 97 ('a') up to 122 ('z')

			Random rand = new Random();
			StringBuilder sb = new StringBuilder(string_length);

			// Build string of random characters
			for (int i = 0; i < string_length; i++)
			{
				byte rand_int = (byte)rand.Next(start, end);    // Generate a random integer
				char monkey_character = (char)rand_int;			// Convert the integer to a lowercase letter
				sb.Append(monkey_character);                    // Add the letter/character to the string object
			}
			return sb.ToString();
		}

		static double CalculateMonkeyProbability(int string_length, int name_length)
		{
			const byte c = 26;                                                  // Number of characters to choose from
			int sequence_combinations = string_length - name_length + 1;        // Number of positions where the name can appear within the random string
			double p = 1 - Math.Pow((1 - Math.Pow(1 / (double)c, name_length))
									 , sequence_combinations);                  // Final probability
			return p;
		}
	}
}