/*
 * Developer's note to Akshay...:
 *		In the Main method, I have put in a total of 5 test cases. The first two are the given tests cases, which I solve via Gaussian Elimination (using all three pivoting strategies).
 * The last three test cases are directly for the Jacobi and Gauss-Seidel methods because I have not put in a function that makes the matrix diagonally dominant first before running
 * the methods, and so the Jacobi and Gauss Seidel fail for the first two test cases and work for the last three. I have also used Gaussian Elimination on the last three test cases
 * to confirm that the Jacobi/Gauss-Seidel solutions are correct.
 */

using System;
using System.Linq;

namespace nvolfango_CA2
{
	class MainClass
	{
		public static void Main()
		{
			// Test cases
			int[,] A1 = new int[,] { { 1, 1, 1, 1 }, { 1, 1, 2, 3 }, { -1, 0, 2, 1 }, { 3, 2, -1, 0 } };
			int[] b1 = new int[] { 1, 2, 1, 1 };

			int[,] A2 = new int[,] { { 3, 1, 4, -1 }, { 2, -2, -1, 2 }, { 5, 7, 14, -8 }, { 1, 3, 2, 4 } };
			int[] b2 = new int[] { 7, 1, 20, -4 };

			GElim system1 = new GElim(A1, b1);
			GElim system2 = new GElim(A2, b2);

			// Test case #1
			Console.WriteLine("Test Case 1:");
			system1.PrintLinearSystem("none");
			system1.PrintLinearSystem("partial");
			system1.PrintLinearSystem("scaled partial");
			Matrix system1_solution = system1.GetSolution("scaled partial");
			Console.WriteLine("\n\nSolution vector for first test case:");
			system1_solution.DisplayMatrix();
			Console.Write("Is there a valid solution? ");
			Console.WriteLine(system1.Solve("scaled partial"));

			Console.WriteLine("\n\n\n");

			// Test case #2
			Console.WriteLine("Test Case 2:");
			system2.PrintLinearSystem("none");
			system2.PrintLinearSystem("partial");
			system2.PrintLinearSystem("scaled partial");
			Matrix system2_solution = system2.GetSolution("scaled partial");
			Console.WriteLine("\n\nSolution vector for first second case:");
			system2_solution.DisplayMatrix();
			Console.Write("Is there a valid solution? ");
			Console.WriteLine(system2.Solve("scaled partial"));

			Console.WriteLine("\n\n\n");

			// Jacobi and Gauss-Seidel Method

			Matrix matrix_1 = new Matrix(new double[,] { { 19, 10, 39 }, { -3, -16, -35 } });
			Matrix matrix_2 = new Matrix(new double[,] { { 10, 2, 3, 23 }, { 5, 25, 7, 76 }, { 9, 10, 32, 125 } });
			Matrix matrix_3 = new Matrix(new double[,] { { 10, 2, 3, 4, 39 }, { 5, 25, 7, 8, 108 }, { 9, 10, 32, 12, 173 }, { 13, 14, 15, 50, 286 } });

			int[,] A3 = new int[,] { { 19, 10 }, { -3, -16 } };
			int[] b3 = new int[] { 39, -35 };

			int[,] A4 = new int[,] { { 10, 2, 3 }, { 5, 25, 7 }, { 9, 10, 32 } };
			int[] b4 = new int[] { 23, 76, 125 };

			int[,] A5 = new int[,] { { 10, 2, 3, 4 }, { 5, 25, 7, 8 }, { 9, 10, 32, 12 }, { 13, 14, 15, 50 } };
			int[] b5 = new int[] { 39, 108, 173, 286 };

			GElim system3 = new GElim(A3, b3);
			GElim system4 = new GElim(A4, b4);
			GElim system5 = new GElim(A5, b5);


			// Test case #3
			Console.WriteLine("\n\nTest Case 3:");
			Matrix Jsolution_matrix1 = Matrix.Solve(matrix_1, "Jacobi");
			Jsolution_matrix1.DisplayMatrix();
			Matrix GSsolution_matrix1 = Matrix.Solve(matrix_1, "Gauss-Seidel");
			GSsolution_matrix1.DisplayMatrix();
			system3.PrintLinearSystem("scaled partial");


			// Test case #4
			Console.WriteLine("\n\nTest Case 4:");
			Matrix Jsolution_matrix2 = Matrix.Solve(matrix_2, "Jacobi");
			Jsolution_matrix2.DisplayMatrix();
			Matrix GSsolution_matrix2 = Matrix.Solve(matrix_2, "Gauss-Seidel");
			GSsolution_matrix2.DisplayMatrix();
			system4.PrintLinearSystem("scaled partial");

			// Test case #5
			Console.WriteLine("\n\nTest Case 5:");
			Matrix Jsolution_matrix3 = Matrix.Solve(matrix_3, "Jacobi");
			Jsolution_matrix3.DisplayMatrix();
			Matrix GSsolution_matrix3 = Matrix.Solve(matrix_3, "Gauss-Seidel");
			GSsolution_matrix3.DisplayMatrix();
			system5.PrintLinearSystem("scaled partial");

		}
	}

	class GElim
	{
		// Properties
		bool solution_exists;
		Matrix solution_vector;					// Vector that will contain the solution matrix/vector for the system of linear equations
		Matrix A;								// Coefficient matrix
		Matrix b;								// Right-hand side vector
		Matrix vector_of_largest_elements;		// Used for G.E. with partial scaling pivoting
		const double zero_tolerance = 1E-12;	// Used to determine if a value is zero
		const double error_tolerance = 1E-10;	// Used to determione if error of a solution is zero

		// Constructors
		public GElim()
		{
			// do nothing
		}

		public GElim(int[,] A, int[] b)
		{
			this.A = new Matrix(A);
			this.b = Matrix.ColumnVector(b);
		}

		// Public Methods
		// Tries to solve the system of linear equations given by A*(solution_vector) = b.
		// Returns true if a satisfactory solution has been found.
		public bool Solve(string pivot_type)
		{
			solution_exists = false;

			if (pivot_type == "scaled partial")
			{
				vector_of_largest_elements = GetLargestElements(A);
			}

			Matrix full = Matrix.Join(A, b);
			int[] pivot_position;
			double pivot_value;
			int num_of_pivots = 0;
			solution_vector = Matrix.ColumnVector(Matrix.ZeroVector(A.ColumnCount));

			// Do Gaussian elimination to reduce the matrix A to upper triangular form.
			for (int c = 0; c < full.ColumnCount; c++)
			{
				pivot_position = FindPivot(ref full, pivot_type, ref num_of_pivots, c, c);

				// If the current column does not have a valid pivot, move to the next column.
				if (c != pivot_position[1])
				{
					continue;
				}

				// If next valid pivot is not in the next row, swap the rows.
				if (num_of_pivots - 1 != pivot_position[0])
				{
					full.SwapRows(num_of_pivots - 1, pivot_position[0]);
					pivot_position[0] = num_of_pivots - 1;
				}

				pivot_value = full.Values[full.RowPos[pivot_position[0]], c];

				// If pivot value is not equal to 1, then divide the row by the pivot_value so that pivot_value = 1.
				if (Math.Abs(pivot_value - 1) > zero_tolerance)
				{
					Matrix.ScaleRow(ref full, full.RowPos[pivot_position[0]], pivot_value);
				}

				// If pivot is not in the last row, check that all values below it in the column are zero.
				// If not, then make them zero.
				if (pivot_position[0] != full.RowCount - 1)
				{
					ReduceBelowPivot(ref full, pivot_position);
				}
			}

			// Do back substitution on the upper triangular matrix to obtain the solution.
			BackSubstitution(ref full, pivot_type);

			// Check if the solution vector is a valid one.
			solution_exists = CheckSolution(A, solution_vector, b);
			Matrix new_solution_vector = new Matrix(solution_vector.RowCount, solution_vector.ColumnCount);

			// Swaps values around in the solution vector in the correct order (due to row swapping).
			for (int r = 0; r < solution_vector.RowCount; r++)
			{
				new_solution_vector.Values[r, 0] = solution_vector.Values[full.RowPos[r], 0];
			}

			solution_vector = new_solution_vector;
			return true;
		}

		public Matrix GetSolution(string pivot_type)
		{
			Solve(pivot_type);
			return solution_vector;
		}

		// Prints the equation consisting of matrices without solving it.
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
						Console.Write("[{0, 3}", A.Values[A.RowPos[r], c]);
					}
					else
					{
						Console.Write("{0, 7}", A.Values[A.RowPos[r], c]);
					}

				}
				Console.Write("] [  x{0}  ]{1,4}{2,5}{3,4}{4,3}", A.RowPos[r], "=", "[", b.Values[A.RowPos[r], 0], "]");
				Console.WriteLine();
			}
		}

		// Prints the equation consisting of matrices, solves it and then prints out the solution.
		public void PrintLinearSystem(string pivot_type)
		{
			Console.WriteLine("\nYour system of equations is as follows:\n");
			for (int r = 0; r < A.RowCount; r++)
			{
				for (int c = 0; c < A.ColumnCount; c++)
				{
					//Console.Write("{0,7}", (Math.Round(values[r, c], 3)));
					if (c == 0)
					{
						Console.Write("[{0, 3}", A.Values[A.RowPos[r], c]);
					}
					else
					{
						Console.Write("{0, 7}", A.Values[A.RowPos[r], c]);
					}

				}
				Console.Write("] [  x{0}  ]{1,4}{2,5}{3,4}{4,3}", A.RowPos[r], "=", "[", b.Values[A.RowPos[r], 0], "]");
				Console.WriteLine();
			}

			solution_exists = Solve(pivot_type);

			if (solution_exists)
			{
				Console.WriteLine("\nWhere your solution matrix (solved via '{0}' pivoting strategy) is given by:\n", pivot_type);
			}
			else
			{
				Console.WriteLine("\n Which does not have a valid solution.\n");
				return;
			}
			
			// Print out solution_vector values
			for (int r = 0; r < A.RowCount; r++)
			{
				Console.Write("[  x{0}  ]", r);
				Console.Write("{0,4}{1,4}{2,3}  ]", "=", "[", solution_vector.Values[r, 0]);
				Console.WriteLine();
			}
		}

		// Private Methods
		// Checks if the solution vector (X) is good enough as a solution to the equation A*X = b
		private bool CheckSolution(Matrix A, Matrix solution_vector, Matrix b)
		{
			double error = Matrix.Norm(Matrix.Subtract(Matrix.Multiply(A, solution_vector), b));
			return (Math.Abs(error) < error_tolerance);
		}

		// Does elementary row operations to turn all values below the pivot to zero.
		private void ReduceBelowPivot(ref Matrix matrix, int[] pivot_position)
		{
			double val;

			for (int r = pivot_position[0] + 1; r < matrix.RowCount; r++)
			{
				val = matrix.Values[matrix.RowPos[r], pivot_position[1]];

				if (Math.Abs(val) > zero_tolerance)
				{
					for (int c = pivot_position[1]; c < matrix.ColumnCount; c++)
					{
						matrix.Values[matrix.RowPos[r], c] -= matrix.Values[matrix.RowPos[pivot_position[0]], c] * val;
					}
				}
			}
		}

		// Finds the next valid pivot in the matrix, starting from the value in the given row and column
		private int[] FindPivot(ref Matrix matrix, string pivot_type, ref int num_of_pivots, int row, int column)
		{
			int[] pivot_position = new int[2];

			if (pivot_type == "none")
			{
				for (int c = column; c < matrix.ColumnCount; c++)
				{
					for (int r = row; r < matrix.RowCount; r++)
					{
						// If current value is zero, continue to the next value.
						if (Math.Abs(matrix.Values[matrix.RowPos[r], column]) < zero_tolerance)
						{
							continue;
						}
						else
						{
							// Once a non-zero value is found, this value will be the next pivot in the elimination.
							num_of_pivots++;
							pivot_position[0] = r;
							pivot_position[1] = c;
							return pivot_position;
						}
					}
				}
			}

			else if (pivot_type == "partial")
			{
				bool pivot_found = false;
				double val, max_val = 0;

				for (int c = column; c < matrix.ColumnCount; c++)
				{
					for (int r = row; r < matrix.RowCount; r++)
					{
						val = matrix.Values[matrix.RowPos[r], column];

						// If current value is zero, continue to the next value.
						if (Math.Abs(val) < zero_tolerance)
						{
							continue;
						}
						else
						{
							// If a non-zero value is found that is the biggest value found so far, the max value is updated to this value
							// and the pivot position is set to be this value.
							if (Math.Abs(val) > Math.Abs(max_val))
							{
								max_val = val;
								pivot_position[0] = r;
								pivot_position[1] = c;
								pivot_found = true;
							}
						}
					}
					// Once all values have been checked, the returned pivot_position will be the greatest value in the column.
					if (pivot_found)
					{
						num_of_pivots++;
						return pivot_position;
					}
				}
			}

			// Note that this is exactly the same as for 'partial', with the only difference being that each value checked in the iteration
			// is first divided by the greatest value along the row of the original coefficient matrix.
			else if (pivot_type == "scaled partial")
			{
				bool pivot_found = false;
				double val, max_val = 0;

				for (int c = column; c < matrix.ColumnCount; c++)
				{
					for (int r = row; r < matrix.RowCount; r++)
					{
						val = matrix.Values[matrix.RowPos[r], column] / vector_of_largest_elements.Values[matrix.RowPos[r], 0];

						// If current value is zero, continue to the next value.
						if (Math.Abs(val) < zero_tolerance)
						{
							continue;
						}
						else
						{
							// If a non-zero value is found that is the biggest value found so far, the max value is updated to this value
							// and the pivot position is set to be this value.
							if (Math.Abs(val) > Math.Abs(max_val))
							{
								max_val = val;
								pivot_position[0] = r;
								pivot_position[1] = c;
								pivot_found = true;
							}
						}
					}
					// Once all values have been checked, the returned pivot_position will be the greatest value in the column.
					if (pivot_found)
					{
						num_of_pivots++;
						return pivot_position;
					}
				}
			}

			pivot_position[0] = pivot_position[1] = -1;
			return pivot_position;
		}


		// Does back subtitution on the given (upper triangular) matrix and populates solution_vector with the correct values.
		private void BackSubstitution(ref Matrix matrix, string pivot_type)
		{
			double solution;
			for (int r = matrix.RowCount - 1; r > -1; r--)
			{
				// The solution is exactly the last value on the same row
				if (r == matrix.RowCount - 1)
				{
					solution_vector.Values[matrix.RowPos[r], 0] = matrix.Values[matrix.RowPos[r], matrix.ColumnCount - 1];
				}
				else
				{
					solution = matrix.Values[matrix.RowPos[r], matrix.ColumnCount - 1];

					for (int c = matrix.ColumnCount - 2; c > r; c--)
					{
						solution -= matrix.Values[matrix.RowPos[r], c] * solution_vector.Values[matrix.RowPos[c], 0];
					}
					solution_vector.Values[matrix.RowPos[r], 0] = solution;
				}
			}
		}

		// This is used when the (Gaussian Elimination) pivoting strategy is 'scaled partial'.
		// This returns a vector of values that are the biggest (in absolute value) values in their respective row in the coefficient matrix.
		private Matrix GetLargestElements(Matrix matrix)
		{
			Matrix vector_of_largest_elements = new Matrix(matrix.RowCount, 1);
			double row_max;

			for (int r = 0; r < matrix.RowCount; r++)
			{
				row_max = 0;
				for (int c = 0; c < matrix.ColumnCount; c++)
				{
					if (matrix.Values[r, c] > row_max)
					{
						row_max = matrix.Values[r, c];
					}
				}
				vector_of_largest_elements.Values[r, 0] = row_max;
			}

			return vector_of_largest_elements;
		}

	}


	class Matrix
	{
		// Properties
		private int nRows, nCols;
		private double[,] values;
		private int[] row_pos;					// Indicates the positions of each row in the matrix. This is used to make row swapping more efficient
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

		// Column-wise joining of matrix X and matrix Y
		public static Matrix Join(Matrix X, Matrix Y)
		{
			Matrix new_matrix = new Matrix(X.nRows, X.nCols + Y.nCols);

			// Fill with values from matrix X
			for (int r = 0; r < X.nRows; r++)
			{
				for (int c = 0; c < X.nCols; c++)
				{
					new_matrix.values[r, c] = X.values[r, c];
				}
			}

			// Fill with values from matrix Y
			for (int r = 0; r < X.nRows; r++)
			{
				for (int c = 0; c < Y.nCols; c++)
				{
					new_matrix.values[r, c + X.nCols] = Y.values[r, c];
				}
			}

			return new_matrix;
		}

		// Prints out values of Matrix object in a nice format
		public void DisplayMatrix()
		{
			Console.WriteLine();
			for (int r = 0; r < nRows; r++)
			{
				for (int c = 0; c < nCols; c++)
				{
					//Console.Write("{0,7}", (Math.Round(values[r, c], 3)));
					Console.Write("{0,7}", values[row_pos[r], c]);
				}
				Console.WriteLine();
			}
			Console.WriteLine();
		}

		// Swaps the row positions of two rows in the Matrix object
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
							Console.WriteLine("number of iterations: {0}", iters);
							Environment.Exit(-1);
						}
					}
					else
					{
						divergence_count = 0;
					}
					x0.values = x1.values;

					iters++;
					error1 = error2;
				}
				while ((iters < max_iters) & (Math.Abs(error2) > error_tolerance));
				Console.WriteLine("\nJacobi solution:");
				Console.WriteLine("number of iterations: {0}", iters);
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
						Console.WriteLine("number of iterations: {0}", iters);
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
				Console.WriteLine("\nGauss-Seidel solution:");
				Console.WriteLine("number of iterations: {0}", iters);
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
	}
}