using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace nvolfango_CA2
{
	class MainClass
	{
		public static void Main(string[] args)
		{
			
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
	}
}