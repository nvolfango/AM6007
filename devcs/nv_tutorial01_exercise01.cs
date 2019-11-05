using System;
using System.Collections.Generic;

namespace nv_tutorial_01
{
	class Program
	{
		public static void Main(string[] args)
		{
			float a, b = 0;
			List<string> operationList = new List<string>(){ "addition", "subtraction", "division", "multiplication" };

			// Ask user for inputs
			Console.WriteLine("Please give me two numbers: ");

			Console.Write("First number: ");
			a = float.Parse(Console.ReadLine());
		
			Console.Write("Second number:");
			b = float.Parse(Console.ReadLine());

			// Ask user for operation they would like to use with the two numbers
			Console.WriteLine("What operation would you like to run on the two numbers?");
			Console.Write("Choices: ");
			operationList.ForEach(operation => Console.WriteLine("{0}\t", operation));
		}
	}
}