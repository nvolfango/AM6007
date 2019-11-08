/* Documentation: Details/Instructions on each method is on the comment above the function definitions
 *  Class: CsvReader
 *		Constructor: CsvReader(string filename, bool header = true, char delimiter = ',')
 *		
 *		Private Methods:
 *			- CsvGetDimensions(string filename, bool header = true, char delimiter = ',')
 *		
 *		Public Methods:
 *			- PrintHeaders()
 *			- PrintValues()
 *			
 *			
 *  Class: LinearEquationSolver
 *		Constructor: LinearEquationSolver(string filename, bool header, char delimiter)
 *		
 *		Private Methods:
 *			- GenerateVariableNames(int variable_count)
 *			- IsInEchelonForm(ref double[] matrix)
 *			- ConvertToEchelonForm(ref double[,] matrix, bool show_steps)
 *			- FindPivotUpper(ref double[,] matrix, int column, ref int pivot_positions_found)
 *			- NormalizeColumnBelowPivot(ref double[,] matrix, int[] pivot_position)
 *			- NormalizeColumnAbovePivot(ref double[,] matrix, int[] pivot_position)
 *			- SwapRows(ref double[,] matrix, int row1, int row2)
 *			- IsConsistent(ref double[,] matrix)
 *			- HasUniqueSolutions(ref double[,] matrix)
 *			
 *		
 *		Public Methods:
 *			- Rref(bool show_steps = false)
 *			- public static void Solve(ref double[,] matrix)
 *			- PrintEquations(ref double[,] matrix)
 *			- PrintMatrix(ref double[,] matrix)
 *			
 *	
 *	Developer's Note:
 *		- Feel free to play around with the Main() function. Make sure to use valid filenames for filename1 and filename2
 */

using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;

namespace CA1
{
	class nvolfango_CA1
	{
		public static void Main(string[] args)
		{

			// Part 1:
			//User inputs:
			string filename1 = @"C:\Users\nvolf\Google Drive 2\5th Year (Masters) Modules\First Semester\AM6007 - Scientific Computing with Numerical Examples - 100% CA\Tutorials\devcs\CA1\test_file_comma.csv";
			bool header1 = true;
			char delimiter1 = ',';

			CsvReader csv_data = new CsvReader(filename1, header1, delimiter1);

			Console.WriteLine("The file location is \"{0}\".", csv_data.Filename);
			Console.WriteLine("The file has {0} columns.", csv_data.ColumnCount);
			Console.WriteLine("The file has {0} rows.", csv_data.RowCount);
			Console.WriteLine("The file has delimiter type \"{0}\".", csv_data.DelimiterType);
			Console.WriteLine("The header option is set to {0}.", csv_data.HasHeader);
			csv_data.PrintHeaders();
			csv_data.PrintValues();
			string[] header_names = csv_data.HeaderNames;


			// Part 2:
			// User inputs:
			string filename2 = @"C:\Users\nvolf\Google Drive 2\5th Year (Masters) Modules\First Semester\AM6007 - Scientific Computing with Numerical Examples - 100% CA\Tutorials\devcs\CA1\test_file_comma.csv";
			bool header2 = false;
			char delimiter2 = ',';
			bool show_steps = false;

			LinearEquationSolver eqs = new LinearEquationSolver(filename2, header2, delimiter2);
			double[,] matrix = eqs.Matrix;
			LinearEquationSolver.PrintEquations(ref matrix);
			double[,] rref_matrix = eqs.Rref(show_steps);
			LinearEquationSolver.PrintEquations(ref rref_matrix);
			LinearEquationSolver.Solve(ref rref_matrix);
		}
	}


	class CsvReader
	{
		/*
		 * Data members
		 */
		string filename;				// Directory and name.extension of the file
		bool header;					// Parameter for specifying inclusion of a header
		string[] header_names;			// Array to store column/header names
		char delimiter;					// Parameter to specify type of delimiter.
		double[,] values;				// Matrix extracted from the file
		int[] dimensions;				// Two-element array (row count, column count)
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
					Console.WriteLine("Warning: {0} values were replaced with a a zero.", num_of_missing_values);
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


	class LinearEquationSolver
	{
		/*
		 * Data members
		 */
		CsvReader csv_data;
		double[,] matrix;				// Matrix of values from file
		int row_count, column_count;
		const double eps = 10E-5;       // Used to find zero values, accounting for
										// numerical inaccuracies resulting from float division.

		/*
		 * Properties
		 */
		public double[,] Matrix
		{
			get
			{
				return matrix;
			}
		}


		/*
		 * Constructor
		 */
		public LinearEquationSolver(string filename, bool header, char delimiter)
		{
			this.csv_data = new CsvReader(filename, header, delimiter);
			this.matrix = this.csv_data.Values;
			this.row_count = csv_data.RowCount;
			this.column_count = csv_data.ColumnCount;
		}


		/*
		 * Private Class Methods - Used by either the class constructor or by one or more of the class public methods
		 */

		// Creates an array of variable names, e.g. 'a', 'b', etc. to be used forvisualizing the system of linear equations.
		// Parameter: variable_count - Number of variables needed to be generated.
		// Note: Since there are 26 letters in the alphabet, the limit of possible variable names is 52 (lowercase + uppercase).
		static char[] GenerateVariableNames(int variable_count)
		{
			int[] lowercase_letters, uppercase_letters;
			char[] variable_names;
			List<int> letters = new List<int>();

			const int length_of_alphabet = 26;
			const int lowercase_start_range = 97; // Integer values that become lowercase letters when casted into char
			const int uppercase_start_range = 65; // Integer values that become uppercase letters when casted into char

			lowercase_letters = Enumerable.Range(lowercase_start_range, length_of_alphabet).ToArray();
			uppercase_letters = Enumerable.Range(uppercase_start_range, length_of_alphabet).ToArray();

			letters.AddRange(lowercase_letters);
			letters.AddRange(uppercase_letters);

			// Convert array of integers into array of letters (lowercase first, then uppercase)
			variable_names = letters.Select(element => (char)element).ToArray();
			variable_names = variable_names.Take(variable_count).ToArray();

			return variable_names;
		}

		// Checks if the matrix is in reduced row echelon form
		bool IsInEchelonForm(double[,] matrix)
		{
			int row_count = matrix.GetLength(0);
			int column_count = matrix.GetLength(1);
			double value;
			bool is_echelon = true; // Assume the matrix is in echelon form, and use for loop to change this if it is not.

			// Variables used to ensure that each successive pivot found is in a different column.
			bool leading_entry_found;
			int col_of_prev_leading_entry;
			int col_of_current_leading_entry = 0;

			// Check each row to make sure that leading entry/pivot is 1 and that all other values in the row are 0.
			for (int i = 0; i < row_count; i++)
			{
				leading_entry_found = false; // Reset this variable to false at the start of each iteration of the row

				if (i == 0)
				{
					col_of_prev_leading_entry = -1;
				}
				else
				{
					col_of_prev_leading_entry = col_of_current_leading_entry;
				}


				for (int j = 0; j < column_count - 1; j++)
				{
					value = matrix[i, j];
					
					// If a value is found that is not zero or one.
					if (Math.Abs(value) > eps && Math.Abs(value - 1) > eps)
					{
						is_echelon = false;
						break;
					}
					else
					{
						if (value == 1)
						{
							// If there is a non-zero value to the right of the leading_entry.
							if (leading_entry_found)
							{
								is_echelon = false;
								break;
							}
							// If this is the first non-zero value found.
							else
							{
								leading_entry_found = true;
								col_of_current_leading_entry = j;
								if (col_of_prev_leading_entry == col_of_current_leading_entry)
								{
									is_echelon = false;
									break;
								}
							}
						}
					}
				}
			}

			return is_echelon;
		}

		// Converts the matrix to echelon form.
		// Parameters: show_steps - If true, prints the system of equations regularly during the conversion process.
		void ConvertToEchelonForm(ref double[,] matrix, bool show_steps)
		{
			bool is_echelon_form;
			int row_count = matrix.GetLength(0);
			int column_count = matrix.GetLength(1);

			is_echelon_form = IsInEchelonForm(matrix);

			int[] pivot_position = new int[2];
			int pivot_positions_found = 0;				// Number of pivots found
			int[] pivot_rows, pivot_columns;			/* Element i of pivot_rows and pivot_columns
														   indicates the position of the ith pivot in the matrix. */

			pivot_rows = pivot_columns = new int[Math.Max(row_count, column_count - 1)];

			if (show_steps) PrintMatrix(ref matrix);

			// Reduce the matrix to upper triangular form
			for (int j = 0; j < column_count - 1; j++)
			{
				if (is_echelon_form) break;

				// Find the next pivot.
				pivot_position = FindPivotUpper(ref matrix, j, ref pivot_positions_found);

				// Normalize pivot row so that pivot value is one, and get the values below the pivot to 0.
				if (pivot_position[0] == pivot_positions_found - 1)
				{
					NormalizeColumnBelowPivot(ref matrix, pivot_position);
				}
				// If next pivot is not in the right position in the matrix, swap some rows.
				else if (pivot_position[0] != -1)
				{
					if (show_steps) Console.WriteLine("Now swapping rows {0} and {1}:", pivot_positions_found, pivot_position[0] + 1);
					SwapRows(ref matrix, pivot_positions_found - 1, pivot_position[0]);
					if (show_steps) PrintMatrix(ref matrix);

					pivot_position[0] = pivot_positions_found - 1; // Accounts for change in pivot position due to swap.
					NormalizeColumnBelowPivot(ref matrix, pivot_position);
				}
				if (show_steps) PrintMatrix(ref matrix);

				pivot_rows[j] = pivot_position[0];
				pivot_columns[j] = pivot_position[1];

				// Regular check to see if matrix reduces to echelon form.
				is_echelon_form = IsInEchelonForm(matrix);
			}

			// Reduce the matrix to a diagonal of ones, as close as possible.
			for (int j = pivot_positions_found - 1; j > -1; j--)
			{
				if (is_echelon_form) break;

				pivot_position[0] = pivot_rows[j];
				pivot_position[1] = pivot_columns[j];

				if (pivot_position[0] != -1)
				{
					NormalizeColumnAbovePivot(ref matrix, pivot_position);
				}

				if (show_steps) PrintMatrix(ref matrix);

				// Regular check to see if matrix reduces to echelon form.
				is_echelon_form = IsInEchelonForm(matrix);
			}
		}

		// Finds the pivot in column of matrix (during the upper triangular matrix conversion process).
		int[] FindPivotUpper(ref double[,] matrix, int column, ref int pivot_positions_found)
		{
			int[] position = new int[2];

			for (int i = column; i < this.row_count; i++)
			{
				if (Math.Abs(matrix[i, column]) > eps)
				{
					position[0] = i;
					position[1] = column;
					pivot_positions_found++;
					return position;
				}
			}
			position[0] = position[0] = -1;
			return position;
		}

		// Changes the values in the pivot row so that pivot value is one, and gets values below the pivot to zero.
		void NormalizeColumnBelowPivot(ref double[,] matrix, int[] pivot_position)
		{
			int pivot_row = pivot_position[0];
			int pivot_column = pivot_position[1];
			double pivot_value = matrix[pivot_row, pivot_column];

			// Scaling pivot row so that pivot point equals 1
			for (int j = pivot_column; j < this.column_count; j++)
			{
				matrix[pivot_row, j] = matrix[pivot_row, j] / pivot_value;
			}

			// Reducing values below the pivot point to 0.
			double scale;
			for (int i = pivot_row + 1; i < this.row_count; i++)
			{
				scale = matrix[i, pivot_column];
				for (int j = pivot_column; j < this.column_count; j++)
				{
					matrix[i, j] -= matrix[pivot_row, j] * scale;
				}
			}
		}

		// Gets values above the pivot to zero.
		void NormalizeColumnAbovePivot(ref double[,] matrix, int[] pivot_position)
		{
			int pivot_row = pivot_position[0];
			int pivot_column = pivot_position[1];
			double pivot_value = matrix[pivot_row, pivot_column];

			// Scaling pivot row so that pivot point equals 1 - bit of overkill, I know...
			for (int j = pivot_column; j > -1; j--)
			{
				matrix[pivot_row, j] = matrix[pivot_row, j] / pivot_value;
			}

			// Reducing values above the pivot point to 0.
			double scale;
			for (int i = pivot_row - 1; i > -1; i--)
			{
				scale = matrix[i, pivot_column];
				for (int j = pivot_column; j < this.column_count; j++)
				{
					matrix[i, j] -= matrix[pivot_row, j] * scale;
				}
			}
		}

		// Switches the values between two rows of a matrix
		void SwapRows(ref double[,] matrix, int row1, int row2)
		{
			double[] tmp = new double[this.column_count];
			for (int j = 0; j < this.column_count; j++)
			{
				tmp[j] = matrix[row1, j];
				matrix[row1, j] = matrix[row2, j];
				matrix[row2, j] = tmp[j];
			}
		}

		// Checks for a row in the matrix where all values are zero except for the last column.
		// If function returns true, then the matrix has at least one solution.
		// If function returns false, then the matrix has no solution.
		static bool IsConsistent(ref double[,] matrix)
		{
			int row_count = matrix.GetLength(0);
			int column_count = matrix.GetLength(1);

			bool row_has_all_zeros = true;

			// Iterate through matrix, starting from bottom-left value
			for (int i = row_count - 1; i > -1; i--)
			{
				for (int j = 0; j < column_count - 1; j++)
				{
					if (Math.Abs(matrix[i, j]) > eps) row_has_all_zeros = false;
				}
				if (row_has_all_zeros && Math.Abs(matrix[i, column_count - 1]) > eps)
				{
					return false;
				}
			}
			return true;
		}

		// Checks if the matrix has a unique set of solution or not.
		static bool HasUniqueSolutions(ref double[,] matrix)
		{
			int row_count = matrix.GetLength(0);
			int column_count = matrix.GetLength(1);
			return (row_count >= column_count - 1);
		}



		// Public Class Methods - Class methods accessible to the user
		public double[,] Rref(bool show_steps = false)
		{
			double[,] solution_matrix = this.matrix;

			PrintMatrix(ref solution_matrix);

			ConvertToEchelonForm(ref solution_matrix, show_steps);

			return solution_matrix;
		}
		
		// Prints out the solution type (no solution, infinite solutions, or unique solution), and
		// prints out a solution if available.
		// Parameters: matrix - Takes the matrix after it has been converted to reduced row echelon format.
		public static void Solve(ref double[,] matrix)
		{
			bool is_consistent;
			bool has_unique_solutions;
			int row_count = matrix.GetLength(0);
			int column_count = matrix.GetLength(1);

			is_consistent = IsConsistent(ref matrix);
			char[] variables_array;

			if (is_consistent)
			{
				has_unique_solutions = HasUniqueSolutions(ref matrix);
				if (has_unique_solutions)
				{
					Console.WriteLine("This particular system of linear equations has a unique solution:");
					variables_array = GenerateVariableNames(column_count - 1);

					for (int j = 0; j < column_count - 1; j++)
					{
						Console.WriteLine("{0} = {1:N4}", variables_array[j], matrix[j, column_count - 1]);
					}

				}
				else
				{
					Console.WriteLine("This particular system of linear equations has infinite solutions.");
					Console.WriteLine("One such sample solution is:");

					variables_array = GenerateVariableNames(column_count - 1);

					double[,] new_matrix = new double[row_count, row_count + 1];
					int num_of_free_variables = (column_count - 1) - row_count;
					double[] free_variables = new double[num_of_free_variables];


					for (int i = 0; i < row_count; i++)
					{
						for (int j = 0; j < num_of_free_variables; j++)
							matrix[i, column_count - 1] -= matrix[i, column_count - j - 2];
					}

					for (int i = 0; i < row_count; i++)
					{
						for (int j = 0; j < row_count; j++)
						{
							new_matrix[i, j] = matrix[i, j];
						}
						new_matrix[i, row_count] = matrix[i, column_count - 1];
					}


					for (int i = 0; i < num_of_free_variables; i++)
					{
						free_variables[i] = 1;
					}

					int k;
					for (k = 0; k < row_count; k++)
					{
						Console.WriteLine("{0} = {1:N4}", variables_array[k], new_matrix[k, row_count]);
					}
					for (int l = 0; l < num_of_free_variables; l++)
					{
						Console.WriteLine("{0} = {1:N4}", variables_array[k + l], free_variables[l]);
					}
					Console.WriteLine("Number of free variables: {0}", num_of_free_variables);
				}
			}
			else
			{
				Console.WriteLine("This particular system of linear equations has no solutions.");
			}

		}

		// Another function to print out the matrix in a nice format, only this is much nicer...
		public static void PrintEquations(ref double[,] matrix)
		{
			int row_count = matrix.GetLength(0);
			int column_count = matrix.GetLength(1);

			const byte max_number_of_variable_names = 52;

			if (column_count - 1 < max_number_of_variable_names)
			{
				char[] variable_names = GenerateVariableNames(column_count - 1);
				for (int i = 0; i < row_count; i++)
				{
					for (int j = 0; j < column_count; j++)
					{
						if (j < column_count - 1)
						{
							double value = matrix[i, j];
							string number;
							char sign;
							char variable_name;

							if (value < 0)
							{
								sign = '-';
								number = value.ToString("#.##").Substring(1);
								variable_name = (char)variable_names[j];
							}
							else if (value > 0)
							{
								sign = '+';
								number = value == 1 ? "\0" : value.ToString("#.##");
								variable_name = (char)variable_names[j];
							}
							else
							{
								sign = '\0';
								number = "\0";
								variable_name = '\0';
							}

							if (j == 0)
							{
								if (Math.Abs(value) < eps)
								{
									Console.Write("{0, 5} ", '\0');
								}
								else
								{
									if (value == 1)
									{
										Console.Write("{0, 5}{1}", '\0', variable_name);
									}
									else
									{
										Console.Write("{0, 5}{1}", value, variable_name);
									}
								}
							}
							else
							{
								Console.Write("{0, 5}{1, 5}{2}", sign, number, variable_name);
							}
						}
						else
						{
							Console.Write("{0, 3}{1, 7:N2}", "=", matrix[i, j]);
						}
					}
					Console.WriteLine();
				}
			}
			else
			{
				Console.WriteLine("Warning: There are too many variables to be printed. Cannot assign a valid (pretty-looking) variable name to each variable.");
			}
			Console.WriteLine();
		}

		// Used to print the values of the matrix in a nice format.
		public static void PrintMatrix(ref double[,] matrix)
		{
			int row_count = matrix.GetLength(0);
			int column_count = matrix.GetLength(1);
			for (int i = 0; i < row_count; i++)
			{
				for (int j = 0; j < column_count; j++)
				{
					Console.Write("{0, 9:N2}", matrix[i, j]);
				}
				Console.WriteLine();
			}
			Console.WriteLine();
		}
	}
}