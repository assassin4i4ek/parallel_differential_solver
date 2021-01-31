package numeric_differential_equation.matrix;

import numeric_differential_equation.solver.Solver;

public interface Matrix {
    double accuracy = 0.0000001;
    int getSize();
    void initSize(int size);
    double get(int i, int j);
    void set(int i, int j, double val);

    default void setRow(double[] row, int i, int offset) {
        for (int j = 0; j < row.length; ++j) {
            set(i,j + offset, row[j]);
        }
    }

    default void parallelSolve(double[] result, int resultOffset, double[] freeTerms, Solver solver) {
        sequentialSolve(result, resultOffset, freeTerms);
    }

    default void sequentialSolve(double[] result, int resultOffset, double[] freeTerms) {
        int size = getSize();
        if (size != freeTerms.length) {
            throw new RuntimeException("Sizes of matrix and free terms don't match");
        }

        for (int k = 0; k < size; ++k) {
            double diagonalElement = get(k, k);
            for (int i = k + 1; i < size; ++i) {
                double iValue = get(i,k);
                if (Math.abs(iValue) > accuracy) {
                    double alpha = iValue / diagonalElement;
                    freeTerms[i] -= alpha * freeTerms[k];
                    for (int j = k;  j < size; ++j) {
                        set(i, j, get(i, j) - alpha * get(k, j));
                    }
                }
            }
        }

        for (int i = size - 1; i >= 0; --i) {
            double sum = 0.0;
            for (int j = i + 1; j < size; ++j) {
                sum += get(i, j) * result[j + resultOffset];
            }
            result[i + resultOffset] = (freeTerms[i] - sum) / get(i, i);
        }
    }

    static Matrix threeDiagonalMatrix() {
        return new ThreeDiagonalMatrix();
    }
}
