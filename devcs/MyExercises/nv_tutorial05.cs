using System;

namespace nv_MatrixClass
{
	class MatrixClass
	{
		// Properties
		private int nRows, nCols;
		private double[,] values;

		// Constructor for matrix of 0's
		public MatrixClass(int nRows, int nCols)
		{
			this.nRows = nRows;
			this.nCols = nCols;
			values = new double[nRows, nCols];
		}

		// Constructor for when the matrix values are known beforehand
		public MatrixClass(double[,] values)
		{
			this.values = values;
			nRows = values.GetLength(0);
			nCols = values.GetLength(1);
		}

		// Getters/Setters
		public int RowCount
		{
			get { return nRows; }
		}

		public int ColumnCount
		{
			get { return nCols; }
		}

		public double[,] Values
		{
			get { return values; }
			set { values = value; }
		}
	

		// Public Methods
		public void DisplayMatrix()
		{
			Console.WriteLine();
			for (int r = 0; r < nRows; r++)
			{
				for (int c = 0; c < nCols; c++)
				{
					Console.Write("{0,5}", values[r, c]);
				}
				Console.WriteLine();
			}
			Console.WriteLine();
		}

		public void Transpose()
		{
			double tmp;

			if (nRows == nCols)
			{
				for (int r = 0; r < nRows; r++)
				{
					for (int c = r + 1; c < nCols; c++)
					{
						tmp = values[r, c];
						values[r, c] = values[c, r];
						values[c, r] = tmp;
					}
				}
			}
			else
			{
				double[,] new_values = new double[nCols, nRows];
				nRows = new_values.GetLength(0);
				nCols = new_values.GetLength(1);

				for (int r = 0; r < nRows; r++)
				{
					for (int c = 0; r < nCols; c++)
					{
						new_values[c, r] = values[r, c];
					}
				}

				values = new_values;
			}
		}

		public void SwapRows(int r1, int r2)
		{
			if (r1 >= nRows)
			{
				throw new ArgumentException("Error: Invalid first argument");
			}
			else if (r2 >= nRows)
			{
				throw new ArgumentException("Error: Invalid second argument");
			}

			double tmp;

			for (int c = 0; c < nCols; c++)
			{
				tmp = values[r1, c];
				values[r1, c] = values[r2, c];
				values[r2, c] = tmp;
			}

		}

		public void SwapColumns(int c1, int c2)
		{
			if (c1 >= nCols)
			{
				throw new ArgumentException("Error: Invalid first argument");
			}
			else if (c2 >= nCols)
			{
				throw new ArgumentException("Error: Invalid second argument");
			}

			double tmp;

			for (int r = 0; r < nRows; r++)
			{
				tmp = values[r, c1];
				values[r, c1] = values[r, c2];
				values[r, c2] = tmp;
			}
		}

		public static MatrixClass Multiply(MatrixClass m1, MatrixClass m2)
		{
			if (m1.nCols != m2.nRows)
			{
				throw new Exception("Error: Matrix dimensions are invalid for matrix multiplication");
			}

			MatrixClass product_matrix = new MatrixClass(m1.nRows, m2.nCols);
			double product_sum;

			for (int r = 0; r < m1.nRows; r++)
			{
				for (int c = 0; c < m2.nCols; c++)
				{
					product_sum = 0;
					for (int i = 0; i < m1.nCols; i++)
					{
						product_sum += m1.Values[r, i] * m2.Values[i, c];
					}
					product_matrix.Values[r, c] = product_sum;
				}
			}

			return product_matrix;
		}

		public static MatrixClass Add(MatrixClass m1, MatrixClass m2)
		{
			if (!((m1.nRows == m2.nRows) & (m1.nCols == m2.nCols)))
			{
				throw new Exception("Matrix dimensions are invalid for matrix addition");
			}

			for (int r = 0; r < m1.nRows; r++)
			{
				for (int c = 0; c < m1.nCols; c++)
				{
					m1.Values[r, c] += m2.Values[r, c];
				}
			}

			return m1;
		}

		public static MatrixClass Subtract(MatrixClass m1, MatrixClass m2)
		{
			if (!((m1.nRows == m2.nRows) & (m1.nCols == m2.nCols)))
			{
				throw new Exception("Matrix dimensions are invalid for matrix addition");
			}

			for (int r = 0; r < m1.nRows; r++)
			{
				for (int c = 0; c < m1.nCols; c++)
				{
					m1.Values[r, c] -= m2.Values[r, c];
				}
			}

			return m1;
		}

		//public static void Main(string[] args)
		//{
		//	//nv_tutorial01_exercise01.MainProgram();
		//	//nv_tutorial01_exercise02.MainProgram();
		//	//nv_tutorial01_exercise03.MainProgram();
		//	//MonteCarloIntegration.MainProgram();

		//	int nRows = 3;
		//	int nCols = 3;

		//	MatrixClass m1 = new MatrixClass(nRows, nCols);
		//	MatrixClass m2 = new MatrixClass(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });

		//	Console.WriteLine("First matrix:");
		//	m1.DisplayMatrix();
		//	Console.WriteLine("Second matrix:");
		//	m2.DisplayMatrix();

		//	// Method tests on m2 matrix

		//	// Row swap
		//	int r1 = 1, r2 = 2;
		//	Console.WriteLine("Matrix before row swap:");
		//	m2.DisplayMatrix();

		//	Console.WriteLine("Swapping row {0} and row {1}", r1, r2);
		//	m2.SwapRows(r1, r2);

		//	Console.WriteLine("Matrix after swap:");
		//	m2.DisplayMatrix();

		//	// Column swap
		//	int c1 = 2, c2 = 0;
		//	Console.WriteLine("Matrix before column swap:");
		//	m2.DisplayMatrix();

		//	Console.WriteLine("Swapping column {0} and column {1}", c1, c2);
		//	m2.SwapColumns(c1, c2);

		//	Console.WriteLine("Matrix after swap:");
		//	m2.DisplayMatrix();

		//	// Transpose
		//	Console.WriteLine("Matrix before transpose:");
		//	m2.DisplayMatrix();

		//	m2.Transpose();

		//	Console.WriteLine("Matrix after transpose:");
		//	m2.DisplayMatrix();

		//	// Multiplication
		//	Console.WriteLine("Matrix before multiplication (with itself):");
		//	m2.DisplayMatrix();

		//	MatrixClass m3 = MatrixClass.Multiply(m2, m2);

		//	Console.WriteLine("Product from matrix multiplication:");
		//	m3.DisplayMatrix();

		//	// Addition
		//	Console.WriteLine("Matrix before addition (with itself):");
		//	m2.DisplayMatrix();

		//	MatrixClass m4 = MatrixClass.Add(m2, m2);

		//	Console.WriteLine("Sum from matrix addition:");
		//	m4.DisplayMatrix();

		//	// Subtraction
		//	Console.WriteLine("Matrix before subtraction (from itself):");
		//	m2.DisplayMatrix();

		//	MatrixClass m5 = MatrixClass.Subtract(m2, m2);

		//	Console.WriteLine("Difference from matrix subtraction:");
		//	m5.DisplayMatrix();
		//}
	}
}

