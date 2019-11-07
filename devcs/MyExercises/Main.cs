using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using nv_tutorial01;
using nv_Topic12_Exercise;
using nv_MatrixClass;

namespace MyExercises
{
	class Program
	{
		public static void Main(string[] args)
		{
			//nv_tutorial01_exercise01.MainProgram();
			//nv_tutorial01_exercise02.MainProgram();
			//nv_tutorial01_exercise03.MainProgram();
			//MonteCarloIntegration.MainProgram();

			int nRows = 3;
			int nCols = 3;

			MatrixClass m1 = new MatrixClass(nRows, nCols);
			MatrixClass m2 = new MatrixClass(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });

			Console.WriteLine("First matrix:");
			m1.DisplayMatrix();
			Console.WriteLine("Second matrix:");
			m2.DisplayMatrix();

			// Method tests on m2 matrix

			// Row swap
			int r1 = 1, r2 = 2;
			Console.WriteLine("Matrix before row swap:");
			m2.DisplayMatrix();

			Console.WriteLine("Swapping row {0} and row {1}", r1, r2);
			m2.SwapRows(r1, r2);

			Console.WriteLine("Matrix after swap:");
			m2.DisplayMatrix();

			// Column swap
			int c1 = 2, c2 = 0;
			Console.WriteLine("Matrix before column swap:");
			m2.DisplayMatrix();

			Console.WriteLine("Swapping column {0} and column {1}", c1, c2);
			m2.SwapColumns(c1, c2);

			Console.WriteLine("Matrix after swap:");
			m2.DisplayMatrix();

			// Transpose
			Console.WriteLine("Matrix before transpose:");
			m2.DisplayMatrix();

			m2.Transpose();

			Console.WriteLine("Matrix after transpose:");
			m2.DisplayMatrix();

			// Multiplication
			Console.WriteLine("Matrix before multiplication (with itself):");
			m2.DisplayMatrix();

			MatrixClass m3 = MatrixClass.Multiply(m2, m2);

			Console.WriteLine("Product from matrix multiplication:");
			m3.DisplayMatrix();

			// Addition
			Console.WriteLine("Matrix before addition (with itself):");
			m2.DisplayMatrix();

			MatrixClass m4 = MatrixClass.Add(m2, m2);

			Console.WriteLine("Sum from matrix addition:");
			m4.DisplayMatrix();

			// Subtraction
			Console.WriteLine("Matrix before subtraction (from itself):");
			m2.DisplayMatrix();

			MatrixClass m5 = MatrixClass.Subtract(m2, m2);

			Console.WriteLine("Difference from matrix subtraction:");
			m5.DisplayMatrix();
		}
	}
}
