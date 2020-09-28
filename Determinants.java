package com.company;

public class Main {


    // testing fitting n points with polynomial of n - 1 th power
    public static void testPolyfit(double [] x, double [] y)
    {
        double [] model = polyfit(x,y);

        double [] results = new double[y.length];

        for(int i = 0; i < x.length; i++)
        {
            double result = 0;
            for(int j = 0; j < model.length; j++)
            {
                result += Math.pow(x[i],j)*model[j];
            }
            results[i] = result;
        }

        for(int i = 0; i< model.length; i++)
        {
            System.out.format("x^(%d) : %f\n",i,model[i]);
        }

        for(int i = 0; i < results.length; i++)
        {
            System.out.format("Fit: %f True: %f\n", results[i],y[i]);
        }
    }


    // fitting n points with n - 1 th power polynomial
    public static double[] polyfit(double[] x, double [] y)
    {
        double [][] matrix = new double[x.length][x.length];

        for(int i = 0; i < x.length; i++)
        {
            for(int j = 0; j < x.length; j++)
            {
                matrix[i][j] = Math.pow(x[i],j);
            }
        }

        double[] coefficients = new double[x.length];

        double mainDeterminant = vandermondeDeterminant(x);

        if(mainDeterminant == 0)
        {
            return coefficients;
        }

        for(int i = 0;i < x.length; i++)
        {
            coefficients[i] = determinant(getSpecificMatrix(matrix,y,i))/mainDeterminant;
        }

        return  coefficients;
    }

    // substituting a column for cramer's rule
    public static double[][] getSpecificMatrix(double [][] matrix, double[] y, int j)
    {
        double [][] result = copyMatrix(matrix);

        for(int i = 0;i< matrix.length;i++)
        {
            result[i][j] = y[i];
        }

        return  result;
    }

    // calculating vandermonder determinant through x vector
    public static double vandermondeDeterminant(double [] x)
    {
        double determinant = 1;
        for(int i = 0; i < x.length; i++)
        {
            for(int j = i + 1; j < x.length; j++)
            {
                determinant = determinant*(x[j]-x[i]);
            }
        }
        return  determinant;
    }

    // recursive function for computing determinants
    public static double recursiveDeterminant(double [][] matrix)
    {
        if(matrix.length==2)
        {
            return matrix[0][0]*matrix[1][1] - matrix[0][1]* matrix[1][0];
        }

        int sign = 1;
        double res = 0;
        for(int i = 0;i < matrix.length;i++)
        {
            res+= sign*matrix[0][i]*recursiveDeterminant(getMinor(matrix,0,i));
            sign = -sign;
        }

        return res;


    }

    // getting a specified minor from a matrix
    public static double[][] getMinor(double[][] matrix,int i, int j)
    {
        double[][] result = new double[matrix.length-1][matrix[0].length-1];
        int iCur =0;
        int jCur = 0;
        for(int k = 0;k< matrix.length;k++)
        {
            for(int m = 0; m< matrix.length;m++)
            {
                if(k!=i)
                {
                    if(m!=j)
                    {
                        result[iCur][jCur] = matrix[k][m];
                        jCur++;
                    }
                }
            }
            jCur = 0;
            if(k!=i)
            {
                iCur++;
            }

        }

        return  result;
    }

    // copying matrix
    public static double[][] copyMatrix(double [][] matrix)
    {
        double[][] copy = new double[matrix.length][];

        for(int i = 0; i<matrix.length; i++)
        {
            copy[i] = new double[matrix[i].length];
            for(int j = 0;j< matrix[i].length;j++)
            {
                copy[i][j] = matrix[i][j];
            }
        }

        return  copy;
    }



    // calculating determinant using gauss's eliminations
    public static double determinant(double[][] a)
    {
        double [][] matrix = copyMatrix(a);
        int n = matrix.length;
        double determinant = 1;

        for(int i = 0;i<n;i++)
        {
            int maxRow = i;

            for(int k = i+1;k<n;k++)
            {
                if(Math.abs(matrix[k][i])>Math.abs(matrix[maxRow][i]))
                {
                    maxRow = k;
                }
            }

            if(matrix[maxRow][i]==0)
            {
                return 0;
            }

            double [] temp = matrix[i];
            matrix[i] = matrix[maxRow];
            matrix[maxRow] = temp;

            determinant*= matrix[i][i];

            if(maxRow!=i)
            {
                determinant = - determinant;
            }

            for(int k = i + 1;k<n;k++)
            {
                double coefficient = matrix[k][i]/matrix[i][i];
                for(int j = i; j<n;j++)
                {
                    matrix[k][j] = matrix[k][j]-matrix[i][j]*coefficient;
                }
            }
        }

        return  determinant;
    }


    public static void main(String[] args)
    {

        System.out.println("" + determinant(new double[][]{{1, 3, 1}, {1, 1, -1}, {3, 11, 5}}));
        System.out.println("" + determinant(new double[][]{{0, 3, 1}, {1, 11, 0}, {3, 0, 5}}));
        System.out.println("" + determinant(new double[][]{{-5, 3, 1}, {1, 10, -9}, {7, -10, 0}}));
        System.out.println("" + determinant(new double[][]{{1, 0, 1}, {5, 1002, 120202}, {0, -8, 1}}));

        System.out.println("" + recursiveDeterminant(new double[][]{{1, 3, 1}, {1, 1, -1}, {3, 11, 5}}));
        System.out.println("" + recursiveDeterminant(new double[][]{{0, 3, 1}, {1, 11, 0}, {3, 0, 5}}));
        System.out.println("" + recursiveDeterminant(new double[][]{{-5, 3, 1}, {1, 10, -9}, {7, -10, 0}}));
        System.out.println("" + recursiveDeterminant(new double[][]{{1, 0, 1}, {5, 1002, 120202}, {0, -8, 1}}));

        double [] x = {1,2,3,4,5};
        double [] y = {-1, 0, 2,5,-1};

        testPolyfit(x,y);

    }
}
