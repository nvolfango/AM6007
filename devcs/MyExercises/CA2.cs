using System;
using System.IO;
using System.Linq;

namespace nvolfango_CA2
{
	class MainClass
	{
		public static void MainProgram()
		{
			bool header = false;
			string ucc_filename = @"Z:\AM6007\My GitHub\devcs\MyExercises\Datasets\test_file_linear_equations5.csv";
			string home_filename = @"C:\Users\nvolf\Google Drive 2\5th Year (Masters) Modules\First Semester\AM6007 - Scientific Computing with Numerical Examples - 100% CA\Tutorials\My GitHub\devcs\MyExercises\Datasets\test_file_linear_equations4.csv";

			CsvReader data = new CsvReader(home_filename, header);
			
			Matrix matrix = new Matrix(data.Values);

			Matrix GEsolution_matrix1 = Matrix.Solve(matrix, "Gaussian Elimination");
			GEsolution_matrix1.DisplayMatrix();

			Matrix GEsolution_matrix2 = Matrix.Solve(matrix, "Gaussian Elimination - partial");
			GEsolution_matrix2.DisplayMatrix();

			Matrix GEsolution_matrix3 = Matrix.Solve(matrix, "Gaussian Elimination - scaled partial");
			GEsolution_matrix3.DisplayMatrix();

			Matrix Jsolution_matrix = Matrix.Solve(matrix, "Jacobi");
			Jsolution_matrix.DisplayMatrix();

			Matrix GSsolution_matrix = Matrix.Solve(matrix, "Gauss-Seidel");
			GSsolution_matrix.DisplayMatrix();		
		}
	}

	class Matrix
	{
		// Properties
		private int nRows, nCols;
		private double[,] values;
		const double zero_tolerance = 1E-16;
		const double error_tolerance = 1E-14;


		// Constructor for matrix of 0's
		public Matrix(int nRows, int nCols)
		{
			this.nRows = nRows;
			this.nCols = nCols;
			values = new double[nRows, nCols];
		}


		// Constructor for when the matrix values are known beforehand
		public Matrix(double[,] values)
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

		// Constructor for identity matrix
		public static Matrix IdentityMatrix(int square_dimension)
		{
			Matrix eye_matrix = new Matrix(square_dimension, square_dimension); ;

			for (int i = 0; i < square_dimension; i++)
			{
				eye_matrix.Values[i, i] = 1;
			}

			return eye_matrix;
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
			if (r1 >= nRows)
			{
				Console.WriteLine("Error: Invalid first argument");
				Environment.Exit(-1);
			}
			else if (r2 >= nRows)
			{
				Console.WriteLine("Error: Invalid second argument");
				Environment.Exit(-1);
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
				Console.WriteLine("Error: Invalid first argument");
				Environment.Exit(-1);
			}
			else if (c2 >= nCols)
			{
				Console.WriteLine("Error: Invalid second argument");
				Environment.Exit(-1);
			}

			double tmp;

			for (int r = 0; r < nRows; r++)
			{
				tmp = values[r, c1];
				values[r, c1] = values[r, c2];
				values[r, c2] = tmp;
			}
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
		public static Matrix Solve(Matrix matrix, string method, int max_iters=2000, int max_divergence_count=5)
		{
			Matrix A = new Matrix(matrix.nRows, matrix.nCols-1);
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
				x0.Values[r, 0] = rand.NextDouble()*20 - 10;
			}

			// Populate matrix b with values from original matrix parameter
			for (int r = 0; r < b.nRows; r++)
			{
				b.Values[r, 0] = matrix.Values[r, matrix.nCols-1];
			}

			if 
			
			else if (method == "Jacobi")
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
			if (matrix.nRows != matrix.nCols-1)
			{
				return false;
			}

			// For each row, calculate absolute value of the diagonal element and of the
			// sum of the other values in the row
			for (int r = 0; r < matrix.nRows; r++)
			{
				absolute_diagonal_value = (ulong) Math.Abs(matrix.Values[r, r]);
				absolute_sum = 0;

				for (int c = 0; c < matrix.nCols-1; c++)
				{
					if (c != r)
					{
						absolute_sum += (ulong) Math.Abs(matrix.Values[r, c]);
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
				for (int c = r + 1; c < matrix.nCols-1; c++)
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
			Matrix inverse_matrix = IdentityMatrix(matrix.nRows);

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

	class CsvReader
	{
		/*
		 * Data members
		 */
		string filename;                // Directory and name.extension of the file
		bool header;                    // Parameter for specifying inclusion of a header
		string[] header_names;          // Array to store column/header names
		char delimiter;                 // Parameter to specify type of delimiter.
		double[,] values;               // Matrix extracted from the file
		int[] dimensions;               // Two-element array (row count, column count)
		int row_count, column_count;


		/*
		 * Properties
		 */

		public double[,] Values
		{
			get
			{
				return values;
			}
		}

		public string[] HeaderNames
		{
			get
			{
				return header_names;
			}
		}

		public bool HasHeader
		{
			get
			{
				return header;
			}
		}

		public string Filename
		{
			get
			{
				return filename;
			}
		}

		public int RowCount
		{
			get
			{
				return row_count;
			}
		}

		public int ColumnCount
		{
			get
			{
				return column_count;
			}
		}

		public char DelimiterType
		{
			get
			{
				return delimiter;
			}
		}


		/*
		 * Constructor
		 */
		public CsvReader(string filename, bool header = true, char delimiter = ',')
		{
			this.filename = filename;
			this.header = header;
			this.delimiter = delimiter;

			// Open file to begin reading data
			Console.WriteLine("Reading file from directory.");
			try
			{
				using (StreamReader sr = new StreamReader(filename))
				{
					string[] line_string;
					double[] line_array;
					int num_of_missing_values = 0;

					// Calculate the dimensions from the file
					dimensions = CsvGetDimensions(filename, header, delimiter);
					row_count = dimensions[0];
					column_count = dimensions[1];

					// Populate the 'header_names' array data member if header is specified true.
					if (header)
					{
						header_names = sr.ReadLine().Split(delimiter);
					}

					// Populate 'values' array with the numerical data from file.
					values = new double[row_count, column_count];
					for (int i = 0; i < row_count; i++)
					{
						// Extract row i from the file into an array (of strings).
						line_string = sr.ReadLine().Split(delimiter);
						num_of_missing_values += (from str in line_string where str == "" select str).Count();

						// Convert to an array of doubles, and replacing missing values with a 0.
						try
						{
							line_array = Array.ConvertAll(line_string, str => str == "" ? 0.0 : double.Parse(str));

							//Fills row i with the values
							for (int j = 0; j < column_count; j++)
							{
								values[i, j] = line_array[j];
							}
						}
						catch (System.FormatException)
						{
							Console.WriteLine("Error: Wrong delimiter used, or there was at least one invalid value encountered.");
							Environment.Exit(-1);
						}
					}
					if (num_of_missing_values > 0)
					{
						Console.WriteLine("Warning: {0} values were replaced with a a zero.", num_of_missing_values);
					}
					Console.WriteLine("The file was read successfully.\n");
				}
			}
			catch (System.IO.FileNotFoundException)
			{
				Console.WriteLine("Error: Specified filename does not exist.");
				Environment.Exit(-1);
			}
		}


		/*
		 * Private Class Methods - Used by either the class constructor or by one or more of the class public methods
		 */

		// Reads through the file specified by filename and calculates the dimensions of the resulting matrix.
		static int[] CsvGetDimensions(string filename, bool header = true, char delimiter = ',')
		{
			int row_count = 0;
			int column_count = 0;
			string[] line;
			int[] dimensions = new int[2];

			using (StreamReader sr_obj = new StreamReader(filename))
			{
				// Iterates through the rows until the end of the file is reached.
				// Column count is the row with the longest count of values
				// Row count is the count of lines read in the file.
				while (!sr_obj.EndOfStream)
				{
					line = sr_obj.ReadLine().Split(delimiter);
					column_count = line.Length > column_count ? line.Length : column_count;
					row_count++;
				}
			}

			// If header = true, the first line of the file is excluded from the row count.
			row_count = header ? row_count - 1 : row_count;

			dimensions[0] = row_count;
			dimensions[1] = column_count;

			return dimensions;
		}


		/*
		 * Public Class Methods - Class methods accessible to the user
		 */

		// Print out the header names from the file, visualised as an array.
		public void PrintHeaders()
		{
			if (header)
			{
				string header_string = String.Join(", ", header_names);
				Console.WriteLine("[" + header_string + "]");
			}
			else
			{
				Console.WriteLine("The file has no headers.");
			}
		}

		// Print out the values from the file in a nice format.
		public void PrintValues()
		{
			// Print out header names if header = true.
			if (header)
			{
				for (int h = 0; h < column_count; h++)
				{
					Console.Write("{0, 6}", header_names[h]);
				}
				Console.WriteLine();
			}
			for (int i = 0; i < row_count; i++)
			{
				for (int j = 0; j < column_count; j++)
				{
					Console.Write("{0, 6}", values[i, j]);
				}
				Console.WriteLine();
			}
		}



	}
}