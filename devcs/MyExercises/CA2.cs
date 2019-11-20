using System;
using System.IO;
using System.Linq;

namespace nvolfango_CA2
{
	class MainClass
	{
		public static void MainProgram()
		{
			int[,] A1 = new int[,] { { 1, 1, 1, 1 }, { 1, 1, 2, 3 }, { -1, 0, 2, 1 }, { 3, 2, -1, 0 } };
			int[] b1 = new int[] { 1, 2, 1, 1 };

			int[,] A2 = new int[,] { { 3, 1, 4, -1 }, { 2, -2, -1, 2 }, { 5, 7, 14, -8 }, { 1, 3, 2, 4 } };
			int[] b2 = new int[] { 7, 1, 20, -4 };

			int[,] A3 = new int[,] { { 3, -1, 3, 1 }, { 6, 0, 9, -2 }, { -12, 0, -10, 5 }, { 72, -8, 48, -19 } };
			int[] b3 = new int[] { 6, 13, -17, 93 };

			GElim system1 = new GElim(A1, b1);
			GElim system2 = new GElim(A2, b2);
			GElim system3 = new GElim(A3, b3);

			//system1.Solve("test");
			system3.Solve("none");
			system3.PrintLinearSystem();
			//system2.PrintLinearSystem();

			//Matrix matrix1 = new Matrix(A1);
			//Matrix matrix1_answer = Matrix.ColumnVector(b1);
			//Matrix matrix2 = new Matrix(A2);
			//Matrix matrix2_answer = Matrix.ColumnVector(b2);
			//matrix1.DisplayMatrix();
			//matrix1_answer.DisplayMatrix();
			//matrix2.DisplayMatrix();
			//matrix2_answer.DisplayMatrix();

			//Matrix GEsolution_matrix1 = Matrix.Solve(matrix, "Gaussian Elimination");
			//GEsolution_matrix1.DisplayMatrix();

			//Matrix GEsolution_matrix2 = Matrix.Solve(matrix, "Gaussian Elimination - partial");
			//GEsolution_matrix2.DisplayMatrix();

			//Matrix GEsolution_matrix3 = Matrix.Solve(matrix, "Gaussian Elimination - scaled partial");
			//GEsolution_matrix3.DisplayMatrix();

			//Matrix Jsolution_matrix = Matrix.Solve(matrix, "Jacobi");
			//Jsolution_matrix.DisplayMatrix();

			//Matrix GSsolution_matrix = Matrix.Solve(matrix, "Gauss-Seidel");
			//GSsolution_matrix.DisplayMatrix();		
		}
	}

	class GElim
	{
		// Properties
		string pivot_type;
		bool solution_exists = false;
		Matrix solution_vector;
		Matrix A;
		Matrix b;
		const double zero_tolerance = 1E-16;
		const double error_tolerance = 1E-14;

		// Constructors
		public GElim()
		{
			// do nothing
		}

		public GElim(int[,] A, int[] b)
		{
			this.A = new Matrix(A);
			this.b = Matrix.ColumnVector(b);
			this.solution_vector = Matrix.ColumnVector(Matrix.ZeroVector(4));
		}

		// Public Methods

		public bool Solve(string pivot_type)
		{
			// Must implement Gaussian elimination followed by back substitution
			Matrix full = Matrix.Join(A, b);
			int[] pivot_position;
			double pivot_value;
			double val;
			int num_of_pivots;

			if (pivot_type == "none")
			{
				for (int c = 0; c < full.ColumnCount; c++)
				{
					num_of_pivots = 0;
					pivot_position = FindPivot(ref full, ref num_of_pivots, c, c);

					// If the current column does not have a valid pivot, move to the next column.
					if (c != pivot_position[1])
					{
						continue;
					}

					for (int r = 0; r < full.RowCount; r++)
					{
						pivot_value = full.Values[full.RowPos[r], c];
						if (Math.Abs(pivot_value - 1) > zero_tolerance)
						{
							Matrix.ScaleRow(ref full, full.RowPos[r], pivot_value);
						}

						// If pivot is not in the last row, check that all values below it in the column are zero.
						// If not, then make them zero.
						if (r != full.RowCount - 1)
						{
							//ReduceBelowPivot(ref full, pivot_position, pivot_value);
							for (int r1 = r + 1; r1 < full.RowCount; r1++)
							{
								val = full.Values[full.RowPos[r1], c];
								if (Math.Abs(val) > 0)
								{
									for (int c1 = c; c1 < full.ColumnCount; c1++)
									{
										full.Values[full.RowPos[r1], c1] -= full.Values[full.RowPos[r1 - 1], c1] * full.Values[full.RowPos[r1], c1];
									}
								}
							}
						}
					}
				}
			}
			else if (pivot_type == "partial")
			{

			}
			else if (pivot_type == "scaled partial")
			{

			}

			return true;
		}

		public Matrix GetSolution()
		{
			return solution_vector;
		}

		public void PrintLinearSystem()
		{
			Console.WriteLine("\nYour system of equations is as follows:\n");
			for (int r = 0; r < A.RowCount; r++)
			{
				for (int c = 0; c < A.ColumnCount; c++)
				{
					//Console.Write("{0,7}", (Math.Round(values[r, c], 3)));
					if (c == 0)
					{
						Console.Write("[{0, 3}", A.Values[r, c]);
					}
					else
					{
						Console.Write("{0, 7}", A.Values[r, c]);
					}

				}
				Console.Write("] [  x{0}  ]{1,4}{2,5}{3,4}{4,3}", r, "=", "[", b.Values[r, 0], "]");
				Console.WriteLine();
			}
			Console.WriteLine("\nWhere your solution matrix is: given by:\n");

			for (int r = 0; r < A.RowCount; r++)
			{
				Console.Write("[  x{0}  ]", r);
				Console.Write("{0,4}{1,4}{2,3}  ]", "=", "[", solution_vector.Values[r, 0]);
				Console.WriteLine();
			}
		}

		public void ReduceBelowPivot(ref Matrix matrix, int[] pivot_position, double pivot_value)
		{
			// If pivot is already in the last row.
			if (pivot_position[0] == matrix.RowCount)
			{
				return;
			}

			for (int r = pivot_position[0]; r < matrix.RowCount; r++)
			{
				for (int c = pivot_position[1]; c < matrix.ColumnCount; c++)
				{

				}
			}
		}

		// Private Methods
		private int[] FindPivot(ref Matrix matrix, ref int num_of_pivots, int row, int column)
		{
			int[] pivot_position = new int[2];
			for (int r = row; r < matrix.RowPos.Length; r++)
			{
				for (int c = column; c < matrix.ColumnCount; c++)
				{
					if (Math.Abs(matrix.Values[matrix.RowPos[r], column]) < zero_tolerance)
					{
						continue;
					}
					else
					{
						num_of_pivots++;
						pivot_position[0] = r;
						pivot_position[1] = c;
						return pivot_position;
					}
				}
			}
			pivot_position[0] = pivot_position[1] = -1;
			return pivot_position;
		}

	}















	class Matrix
	{
		// Properties
		private int nRows, nCols;
		private double[,] values;
		private int[] row_pos;
		const double zero_tolerance = 1E-16;
		const double error_tolerance = 1E-14;


		// Getters/Setters
		public int RowCount
		{
			get { return nRows; }
		}

		public int ColumnCount
		{
			get { return nCols; }
		}

		public int[] RowPos
		{
			get { return row_pos; }
		}

		public double[,] Values
		{
			get { return values; }
			set { values = value; }
		}


		// Constructor with no inputs
		public Matrix()
		{
			// do nothing
		}

		// Constructor for matrix of 0's
		public Matrix(int nRows, int nCols)
		{
			this.nRows = nRows;
			this.nCols = nCols;
			values = new double[nRows, nCols];
			row_pos = Enumerable.Range(0, nRows).ToArray();
		}


		// Constructor for when the matrix values are known beforehand (given as array of doubles)
		public Matrix(double[,] values)
		{
			this.values = values;
			nRows = values.GetLength(0);
			nCols = values.GetLength(1);
			row_pos = Enumerable.Range(0, nRows).ToArray();
		}

		// Constructor for when the matrix values are known beforehand (given as array of ints)
		public Matrix(int[,] values)
		{
			double[,] double_values = new double[values.GetLength(0), values.GetLength(1)];
			for (int r = 0; r < values.GetLength(0); r++)
			{
				for (int c = 0; c < values.GetLength(1); c++)
				{
					double_values[r, c] = values[r, c];
				}
			}
			this.values = double_values;
			nRows = values.GetLength(0);
			nCols = values.GetLength(1);
			row_pos = Enumerable.Range(0, nRows).ToArray();
		}


		
		// Public Methods

		// Zero vector
		public static double[] ZeroVector(int length)
		{
			return new double[length];
		}

		// Identity matrix
		public static Matrix IdentityMatrix(int square_dimension)
		{
			Matrix eye_matrix = new Matrix(square_dimension, square_dimension);

			for (int i = 0; i < square_dimension; i++)
			{
				eye_matrix.Values[i, i] = 1;
			}

			return eye_matrix;
		}

		// Row vector
		public static Matrix RowVector(double[] values)
		{
			Matrix row_vector = new Matrix(1, values.Length);
			
			for (int c = 0; c < values.Length; c++)
			{
				row_vector.values[0, c] = values[c];
			}

			return row_vector;
		}

		public static Matrix RowVector(int[] values)
		{
			Matrix row_vector = new Matrix(1, values.Length);

			for (int c = 0; c < values.Length; c++)
			{
				row_vector.values[0, c] = values[c];
			}

			return row_vector;
		}

		// Column vector
		public static Matrix ColumnVector(double[] values)
		{
			Matrix column_vector = new Matrix(values.Length, 1);

			for (int r = 0; r < values.Length; r++)
			{
				column_vector.values[r, 0] = values[r];
			}

			return column_vector;
		}

		public static Matrix ColumnVector(int[] values)
		{
			Matrix column_vector = new Matrix(values.Length, 1);

			for (int r = 0; r < values.Length; r++)
			{
				column_vector.values[r, 0] = values[r];
			}

			return column_vector;
		}


		public static Matrix Join(Matrix X, Matrix Y)
		{
			Matrix new_matrix = new Matrix(X.nRows, X.nCols + Y.nCols);

			for (int r = 0; r < X.nRows; r++)
			{
				for (int c = 0; c < X.nCols; c++)
				{
					new_matrix.values[r, c] = X.values[r, c];
				}
			}

			for (int r = 0; r < X.nRows; r++)
			{
				for (int c = 0; c < Y.nCols; c++)
				{
					new_matrix.values[r, c + X.nCols] = Y.values[r, c];
				}
			}

			return new_matrix;
		}


		public void DisplayMatrix()
		{
			Console.WriteLine();
			for (int r = 0; r < nRows; r++)
			{
				for (int c = 0; c < nCols; c++)
				{
					//Console.Write("{0,7}", (Math.Round(values[r, c], 3)));
					Console.Write("{0,7}", values[r, c]);
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
			int temp = row_pos[r1];
			row_pos[r1] = row_pos[r2];
			row_pos[r2] = temp;
		}


		public static Matrix Multiply(Matrix m1, Matrix m2)
		{
			if (m1.nCols != m2.nRows)
			{
				Console.WriteLine("Error: Matrix dimensions are invalid for matrix multiplication");
				Environment.Exit(-1);
			}

			Matrix product_matrix = new Matrix(m1.nRows, m2.nCols);
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


		public static Matrix Add(Matrix m1, Matrix m2)
		{
			if (!((m1.nRows == m2.nRows) & (m1.nCols == m2.nCols)))
			{
				Console.WriteLine("Matrix dimensions are invalid for matrix addition");
				Environment.Exit(-1);
			}

			Matrix sum_matrix = new Matrix(m1.nRows, m1.nCols);

			for (int r = 0; r < m1.nRows; r++)
			{
				for (int c = 0; c < m1.nCols; c++)
				{
					sum_matrix.Values[r, c] = m1.Values[r, c] + m2.Values[r, c];
				}
			}

			return sum_matrix;
		}


		public static Matrix Subtract(Matrix m1, Matrix m2)
		{
			if (!((m1.nRows == m2.nRows) & (m1.nCols == m2.nCols)))
			{
				Console.WriteLine("Matrix dimensions are invalid for matrix subtraction");
				Environment.Exit(-1);
			}

			Matrix difference_matrix = new Matrix(m1.nRows, m1.nCols);

			for (int r = 0; r < m1.nRows; r++)
			{
				for (int c = 0; c < m1.nCols; c++)
				{
					difference_matrix.Values[r, c] = m1.Values[r, c] - m2.Values[r, c];
				}
			}

			return difference_matrix;
		}


		// The assumption is that the matrix parameter is given as an m x (n+1) matrix, where each row
		// is a linear equation in n variables, with the last column being the answer.
		// The equation, in matrix form, is given by: Ax = b, where
		//		A is a matrix consisting of the first n columns of the matrix parameter
		//		x is the vector of variables
		//		b is the last column of the matrix parameter
		public static Matrix Solve(Matrix matrix, string method, int max_iters = 2000, int max_divergence_count = 5)
		{
			Matrix A = new Matrix(matrix.nRows, matrix.nCols - 1);
			Matrix x0 = new Matrix(matrix.nRows, 1);
			Matrix x1 = new Matrix(matrix.nRows, 1);
			Matrix b = new Matrix(matrix.nRows, 1);
			int iters = 0;
			int divergence_count = 0;

			// Populate matrix A with values from original matrix parameter
			for (int r = 0; r < A.nRows; r++)
			{
				for (int c = 0; c < A.nCols; c++)
				{
					A.Values[r, c] = matrix.Values[r, c];
				}
			}

			// Assign random initial values to matrix x 
			Random rand = new Random();
			for (int r = 0; r < x0.nRows; r++)
			{
				// Random initial values between -10 and 10
				x0.Values[r, 0] = rand.NextDouble() * 20 - 10;
			}

			// Populate matrix b with values from original matrix parameter
			for (int r = 0; r < b.nRows; r++)
			{
				b.Values[r, 0] = matrix.Values[r, matrix.nCols - 1];
			}

			//if 

			if (method == "Jacobi")
			{
				Matrix R = OffDiagonals(A);
				Matrix D = Diagonals(A);
				Matrix D_inv = InvertDiagonalMatrix(D);
				Matrix D_inv_R = Multiply(D_inv, R);

				double error1 = 0;
				double error2;

				do
				{
					x1 = Multiply(D_inv, Subtract(b, Multiply(R, x0)));
					error2 = Norm(Subtract(Multiply(A, x1), b));
					if (error2 > error1)
					{
						divergence_count++;
						if (divergence_count == max_divergence_count)
						{
							Console.WriteLine("Error: The solution is diverging.");
							Environment.Exit(-1);
						}
					}
					x0.values = x1.values;

					iters++;
					error1 = error2;
				}
				while ((iters < max_iters) & (Math.Abs(error2) > error_tolerance));
				Console.WriteLine("iters = {0}", iters);
				Console.WriteLine("error: {0}", error2);
			}
			else if (method == "Gauss-Seidel")
			{
				Matrix L = Triangle(A, "Lower");
				Matrix D = Diagonals(A);
				Matrix U = Triangle(A, "Upper");

				Matrix Ld = Add(L, D);
				Matrix Ld_inv = InvertLowerTriangularMatrix(Ld);

				double error1 = 0;
				double error2;

				do
				{
					x1 = Multiply(Ld_inv, Subtract(b, Multiply(U, x0)));
					error2 = Norm(Subtract(Multiply(A, x1), b));
					if (error2 > error1)
					{
						divergence_count++;
					}
					else
					{
						divergence_count = 0;
					}
					if (divergence_count == max_divergence_count)
					{
						Console.WriteLine("Error: The solution is diverging.");
						Console.WriteLine("iters: {0}", iters);
						Environment.Exit(-1);
					}
					x0.Values = x1.Values;
					iters++;
					error1 = error2;

					if (iters == max_iters - 1)
					{
						Console.WriteLine("Warning: maximum number of iterations ({0}) reached.", max_iters);
					}
				}
				while (iters < max_iters & Math.Abs(error2) > zero_tolerance);
				Console.WriteLine("iters = {0}", iters);
				Console.WriteLine("error: {0}", error2);
			}

			return x1;
		}


		// Returns a matrix consisting of the diagonals of the input matrix,
		// with off-diagonal elements set to zero.
		public static Matrix Diagonals(Matrix matrix)
		{
			Matrix diagonal_matrix = new Matrix(matrix.nRows, matrix.nCols);
			int min_dimension = Math.Min(matrix.nRows, matrix.nCols);

			for (int i = 0; i < min_dimension; i++)
			{
				diagonal_matrix.Values[i, i] = matrix.Values[i, i];
			}

			return diagonal_matrix;
		}


		// Returns a matrix consisting of the the upper/lower triangle
		// of the input matrix, with other elements set to zero.
		public static Matrix Triangle(Matrix matrix, string shape)
		{
			Matrix triangle = new Matrix(matrix.nRows, matrix.nCols);

			if (shape == "Lower")
			{
				for (int r = 1; r < matrix.nRows; r++)
				{
					for (int c = 0; c < r; c++)
					{
						triangle.Values[r, c] = matrix.Values[r, c];
					}
				}
			}
			else if (shape == "Upper")
			{
				for (int c = 1; c < matrix.nCols; c++)
				{
					for (int r = 0; r < c; r++)
					{
						triangle.Values[r, c] = matrix.Values[r, c];
					}
				}
			}

			return triangle;
		}


		// Returns a matrix consisting of the off-diagonal vlues of the input matrix,
		// with elements on the diagonal set to zero.
		public static Matrix OffDiagonals(Matrix matrix)
		{
			Matrix off_diagonal_matrix = new Matrix(matrix.nRows, matrix.nCols);

			for (int r = 0; r < matrix.nRows; r++)
			{
				for (int c = 0; c < matrix.nCols; c++)
				{
					if (r != c)
					{
						off_diagonal_matrix.Values[r, c] = matrix.Values[r, c];
					}
				}
			}

			return off_diagonal_matrix;
		}


		// Returns a matrix consisting of the inverse of a diagonal matrix.
		public static Matrix InvertDiagonalMatrix(Matrix matrix)
		{
			Matrix inverse_matrix = new Matrix(matrix.nRows, matrix.nCols);

			for (int i = 0; i < matrix.nCols; i++)
			{
				if (Math.Abs(matrix.Values[i, i]) < zero_tolerance)
				{
					Console.WriteLine("Error: Matrix is not invertible. Cannot find a solution via the Jacobi method.");
					Environment.Exit(-1);
				}
			}

			for (int i = 0; i < matrix.nRows; i++)
			{
				inverse_matrix.Values[i, i] = 1 / matrix.Values[i, i];
			}

			return inverse_matrix;
		}


		public static bool IsDiagonallyDominant(Matrix matrix)
		{
			ulong absolute_diagonal_value;
			ulong absolute_sum;

			// Check if matrix is square
			if (matrix.nRows != matrix.nCols - 1)
			{
				return false;
			}

			// For each row, calculate absolute value of the diagonal element and of the
			// sum of the other values in the row
			for (int r = 0; r < matrix.nRows; r++)
			{
				absolute_diagonal_value = (ulong)Math.Abs(matrix.Values[r, r]);
				absolute_sum = 0;

				for (int c = 0; c < matrix.nCols - 1; c++)
				{
					if (c != r)
					{
						absolute_sum += (ulong)Math.Abs(matrix.Values[r, c]);
					}
				}

				// If a row is found that fails the criterion, then matrix is not diagonally dominant.
				if (absolute_diagonal_value < absolute_sum)
				{
					return false;
				}
			}

			// If the algorithm reaches this point, then it remains that the matrix is diagonally dominant.
			return true;
		}


		//Calculates the norm of a vector
		public static double Norm(Matrix vector)
		{
			double norm;
			double sum = 0;

			for (int r = 0; r < vector.nRows; r++)
			{
				sum += vector.Values[r, 0] * vector.Values[r, 0];
			}

			norm = Math.Sqrt(sum);

			return norm;
		}


		public static bool IsSymmetric(Matrix matrix)
		{
			// Check if matrix is square
			if (matrix.nRows != matrix.nCols - 1)
			{
				return false;
			}

			// For each value in the upper triangle, check if it is equal to 
			//the corresponding value in the lower triangle.
			for (int r = 0; r < matrix.nRows; r++)
			{
				for (int c = r + 1; c < matrix.nCols - 1; c++)
				{
					if (matrix.Values[r, c] != matrix.Values[c, r])
					{
						return false;
					}
				}
			}

			// If no 'unequal' terms are found, the obvious result is that the matrix is symmetric.
			return true;
		}


		// Calculates the inverse of a matrix, assuming the matrix is already lower triangular form.
		// Returns the calculated inverse matrix.
		public static Matrix InvertLowerTriangularMatrix(Matrix matrix)
		{
			Matrix inverse_matrix = Matrix.IdentityMatrix(matrix.nRows);

			for (int i = 0; i < matrix.nCols; i++)
			{
				if (Math.Abs(matrix.Values[i, i]) < zero_tolerance)
				{
					Console.WriteLine("Error: Matrix is not invertible. Cannot find a solution via the Gauss-Seidel method.");
					Environment.Exit(-1);
				}
			}

			for (int c = 0; c < matrix.nCols; c++)
			{
				double diagonal_value = matrix.Values[c, c];

				if (diagonal_value != 1)
				{
					ScaleRow(ref inverse_matrix, c, diagonal_value);
				}

				for (int r = c + 1; r < matrix.nRows; r++)
				{
					for (int c1 = 0; c1 < c + 1; c1++)
					{
						inverse_matrix.Values[r, c1] -= matrix.Values[r, c] * inverse_matrix.Values[c, c1];
					}
				}
			}

			return inverse_matrix;
		}


		public static void ScaleRow(ref Matrix matrix, int row, double scale_value)
		{
			for (int c = 0; c < matrix.nCols; c++)
			{
				matrix.Values[row, c] /= scale_value;
			}
		}

		//public static Matrix Submatrix(Matrix matrix, int excluded_row, int excluded_column)
		//{
		//	Matrix submatrix = new Matrix(matrix.nRows-1, matrix.nCols-1);
		//	int r1 = 0;
		//	int c1;

		//	for (int r = 0; r < matrix.nRows; r++)
		//	{
		//		c1 = 0;
		//		for (int c = 0; c < matrix.nCols; c++)
		//		{
		//			if ((r == excluded_row) | (c == excluded_column))
		//			{
		//				continue;
		//			}
		//			else
		//			{
		//				submatrix.Values[r1, c1] = matrix.Values[r, c];
		//			}
		//			if (c != excluded_column) c1++;
		//		}
		//		if (r != excluded_row) r1++;
		//	}

		//	return submatrix;
		//}

		//public static double Determinant(Matrix matrix)
		//{
		//	if (matrix.nRows != matrix.nCols)
		//	{
		//		Console.WriteLine("Error: Matrix has no determinant.");
		//		Environment.Exit(-1);
		//	}
		//	double determinant = 0;
		//	double coefficient;

		//	if *matrix.nRows 

		//	for (int r = 0; r < matrix.nRows; r++)
		//	{
		//		for (int c = 0; c < matrix.nCols; c++)
		//		{
		//			coefficient = c % 2 == 0 ? matrix.Values[r, c] : -matrix.Values[r, c];
		//			return (coefficient * Determinant(Submatrix(matrix, r, c)));
		//		}
		//	}

		//	return determinant;
		//}
	}
}