using System;
using System.IO;
using System.Linq;

namespace nvolfango_CA2
{
	class MainClass
	{
		public static void MainProgram()
		{
			CsvReader data = new CsvReader(@"Z:\AM6007\My GitHub\devcs\MyExercises\Datasets\test_file_linear_equations.csv", header:false);
			
			Matrix matrix = new Matrix(data.Values);
			matrix.DisplayMatrix();
			
		}
	}

	class Matrix
	{
		// Properties
		private int nRows, nCols;
		private double[,] values;


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

		public static Matrix Multiply(Matrix m1, Matrix m2)
		{
			if (m1.nCols != m2.nRows)
			{
				throw new Exception("Error: Matrix dimensions are invalid for matrix multiplication");
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

		public static Matrix Subtract(Matrix m1, Matrix m2)
		{
			if (!((m1.nRows == m2.nRows) & (m1.nCols == m2.nCols)))
			{
				throw new Exception("Matrix dimensions are invalid for matrix subtraction");
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


		// The assumption is that the matrix parameter is given as an m x (n+1) matrix, where each row
		// is a linear equation in n variables, with the last column being the answer.
		// The equation, in matrix form, is given by: Ax = b, where
		//		A is a matrix consisting of the first n columns of the matrix parameter
		//		x is the vector of variables
		//		b is the last column of the matrix parameter
		public Matrix Solve(string method, Matrix matrix, double tolerance=1E-2, int max_iters=1000)
		{
			Matrix A = new Matrix(matrix.nRows, matrix.nCols-1);
			Matrix x0 = new Matrix(matrix.nRows, 1);
			Matrix x1 = new Matrix(matrix.nRows, 1);
			Matrix b = new Matrix(matrix.nRows, 1);
			int iters = 0;

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
				x0.Values[r, 1] = rand.NextDouble()*20 - 10;
			}

			// Populate matrix b with values from original matrix parameter
			for (int r = 0; r < b.nRows; r++)
			{
				b.Values[r, 1] = matrix.Values[r, matrix.nCols-1];
			}


			if (method == "Jacobi")
			{
				Matrix R = OffDiagonals(A);
				Matrix D = Diagonals(A);
				Matrix D_inv = InvertDiagonalMatrix(D);
				double error;

				do
				{
					x1 = Multiply(D_inv, Subtract(b, Multiply(R, x0)));
					error = Norm(Subtract(Multiply(A, x1), b));
					x0 = x1;
					iters++;
				}
				while (iters < max_iters & error < tolerance);

			}
			else if (method == "Gauss-Seidel")
			{
				Matrix R = OffDiagonals(A);
				Matrix D = Diagonals(A);
				
			}

			return x1;
		}


		// Private Methods

		// Returns a matrix consisting of the diagonals of the input matrix,
		// with off-diagonal elements set to zero.
		Matrix Diagonals(Matrix matrix)
		{
			Matrix diagonal_matrix = new Matrix(matrix.nRows, matrix.nCols);
			int min_dimension = Math.Min(matrix.nRows, matrix.nCols);

			for (int i = 0; i < min_dimension; i++)
			{
				diagonal_matrix.Values[i, i] = matrix.Values[i, i];
			}

			return diagonal_matrix;
		}

		// Returns a matrix consisting of the off-diagonal vlues of the input matrix,
		// with elements on the diagonal set to zero.
		Matrix OffDiagonals(Matrix matrix)
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
		Matrix InvertDiagonalMatrix(Matrix matrix)
		{
			int min_dimension = Math.Min(matrix.nRows, matrix.nCols);
			for (int i = 0; i < min_dimension; i++)
			{
				matrix.Values[i, i] = 1 / matrix.Values[i, i];
			}

			return matrix;
		}

		bool IsDiagonallyDominant(Matrix matrix)
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
		double Norm(Matrix vector)
		{
			double norm;
			double sum = 0;

			for (int r = 0; r < vector.nRows; r++)
			{
				sum += vector.Values[r, 1] * vector.Values[r, 1];
			}

			norm = Math.Sqrt(sum);
			
			return norm;
		}

		//bool IsPostiveDefinite(Matrix matrix)
		//{


		//	return true;
		//}

		bool IsSymmetric(Matrix matrix)
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

		Matrix Rref(Matrix matrix)
		{

		}

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