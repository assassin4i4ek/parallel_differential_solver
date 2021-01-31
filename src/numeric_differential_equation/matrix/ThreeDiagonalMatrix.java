package numeric_differential_equation.matrix;

import mpi.MPI;
import numeric_differential_equation.solver.Solver;

import java.util.ArrayList;
import java.util.List;

import static numeric_differential_equation.mpi_features.MpiTags.*;

public class ThreeDiagonalMatrix implements Matrix {
    private int size = 0;
    private List<Double> mainDiagonal = new ArrayList<>();
    private List<Double> upperDiagonal = new ArrayList<>();
    private List<Double> lowerDiagonal = new ArrayList<>();
    private double additionalLeftElement = 0.0;
    private double additionalRightElement = 0.0;

    @Override
    public int getSize() {
        return size;
    }

    @Override
    public void initSize(int size) {
        for (int i = mainDiagonal.size(); i < size; ++i) {
            mainDiagonal.add(0.0);
        }
        for (int i = upperDiagonal.size(); i < size - 1; ++i) {
            upperDiagonal.add(0.0);
        }
        for (int i = lowerDiagonal.size(); i < size - 1; ++i) {
            lowerDiagonal.add(0.0);
        }
        this.size = size;
    }

    @Override
    public double get(int i, int j) {
        if (i == j) {
            return mainDiagonal.get(i);
        }
        else if (i == j - 1) {
            return upperDiagonal.get(i);
        }
        else if (i == j + 1) {
            return lowerDiagonal.get(j);
        }
        else {
            return 0.0;
        }
    }

    @Override
    public void set(int i, int j, double val) {
        if (i == 0 && j == -1) {
            additionalLeftElement = val;
        }
        else if (i == size - 1 && j == size) {
            additionalRightElement = val;
        }
        else {
            if (i == j) {
                for (int count = mainDiagonal.size(); count <= i; ++count) {
                    mainDiagonal.add(0.0);
                }
                mainDiagonal.set(i, val);
            }
            else if (i == j - 1) {
                for (int count = upperDiagonal.size(); count <= i; ++count) {
                    upperDiagonal.add(0.0);
                }
                upperDiagonal.set(i, val);
            }
            else if (i == j + 1) {
                for (int count = lowerDiagonal.size(); count <= j; ++count) {
                    lowerDiagonal.add(0.0);
                }
                lowerDiagonal.set(j, val);
            }
            else {
                throw new RuntimeException("Matrix must be three diagonal. Indexes: " + i + ", " + j);
            }
            if (i + 1 > size) {
                size = i + 1;
                if (additionalRightElement != 0.0) {
                    upperDiagonal.add(additionalRightElement);
                    additionalRightElement = Double.NaN;
                }
            }
            if (j + 1 > size) {
                size = j + 1;
            }
        }
    }

    @Override
    public void parallelSolve(double[] result, int resultOffset, double[] freeTerms, Solver solver) {
        int rank = solver.getRank();
        int numOfProcesses = solver.getNumOfProcesses();
        int count = solver.getCount(rank);
        double leftColumn[] = new double[0];
        double rightColumn[] = new double[0];
        if (rank != 0) {
            leftColumn = new double[size];
            leftColumn[0] = additionalLeftElement;
        }
        if (rank != numOfProcesses - 1) {
            rightColumn = new double[size];
            rightColumn[rightColumn.length - 1] = additionalRightElement;
        }

        for (int i = 1; i < size; ++i) {
            double alpha = lowerDiagonal.get(i - 1) / mainDiagonal.get(i - 1);
            lowerDiagonal.set(i - 1, 0.0);
            mainDiagonal.set(i, mainDiagonal.get(i) - alpha * upperDiagonal.get(i - 1));
            freeTerms[i] -= alpha * freeTerms[i - 1];
            if (rank != 0) {
                leftColumn[i] = - alpha * leftColumn[i - 1];
            }
        }

        for (int i = size - 2; i >= 0; --i) {
            double beta = upperDiagonal.get(i) / mainDiagonal.get(i + 1);
            upperDiagonal.set(i, 0.0);
            freeTerms[i] -= beta * freeTerms[i + 1];
            if (rank != 0) {
                leftColumn[i] -= beta * leftColumn[i + 1];
            }
            if (rank != numOfProcesses - 1) {
                rightColumn[i] = - beta * rightColumn[i + 1];
            }
        }

        double buf[] = new double[1];
        double upperBoundResults[] = new double[rank == 0 ? 1 : 2];
        double lowerBoundResults[] = new double[rank == numOfProcesses - 1 ? 1 : 2];

        if (rank == 0) {
            ThreeDiagonalMatrix newThreeDiagonalMatrix = new ThreeDiagonalMatrix();
            List<Double> newFreeTerms = new ArrayList<>();

            newThreeDiagonalMatrix.mainDiagonal.add(mainDiagonal.get(0));
            newThreeDiagonalMatrix.upperDiagonal.add(rightColumn[0]);
            newFreeTerms.add(freeTerms[0]);

            if (count > 1) {
                newThreeDiagonalMatrix.upperDiagonal.add(mainDiagonal.get(size - 1));
                newThreeDiagonalMatrix.mainDiagonal.add(rightColumn[size - 1]);
                newThreeDiagonalMatrix.lowerDiagonal.add(0.0);
                newFreeTerms.add(freeTerms[size - 1]);
            }

            for (int i = 1; i < numOfProcesses; ++i) {
                int iCount = solver.getCount(i) - (i == numOfProcesses - 1 ? 1 : 0);
                if (iCount > 1) {
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, MAIN_DIAGONAL_FIRST_ELEMENT);
                    newThreeDiagonalMatrix.lowerDiagonal.add(buf[0]);
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, LEFT_COLUMN_FIRST_ELEMENT);
                    newThreeDiagonalMatrix.mainDiagonal.add(buf[0]);
                    if (i != numOfProcesses - 1) {
                        MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, RIGHT_COLUMN_FIRST_ELEMENT);
                        newThreeDiagonalMatrix.upperDiagonal.add(buf[0]);
                    }
                    else {
                        newThreeDiagonalMatrix.upperDiagonal.add(0.0);
                    }
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, LEFT_COLUMN_LAST_ELEMENT);
                    newThreeDiagonalMatrix.lowerDiagonal.add(buf[0]);
                    if (i != numOfProcesses - 1) {
                        MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, RIGHT_COLUMN_LAST_ELEMENT);
                        newThreeDiagonalMatrix.mainDiagonal.add(buf[0]);
                        MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, MAIN_DIAGONAL_LAST_ELEMENT);
                        newThreeDiagonalMatrix.upperDiagonal.add(buf[0]);
                    }
                    else {
                        MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, MAIN_DIAGONAL_LAST_ELEMENT);
                        newThreeDiagonalMatrix.mainDiagonal.add(buf[0]);
                    }
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, FIRST_FREE_TERM);
                    newFreeTerms.add(buf[0]);
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, LAST_FREE_TERM);
                    newFreeTerms.add(buf[0]);
                }
                else if (iCount == 1) {
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, MAIN_DIAGONAL_FIRST_ELEMENT);
                    newThreeDiagonalMatrix.lowerDiagonal.add(buf[0]);
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, LEFT_COLUMN_FIRST_ELEMENT);
                    newThreeDiagonalMatrix.mainDiagonal.add(buf[0]);
                    if (i != numOfProcesses - 1) {
                        MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, RIGHT_COLUMN_FIRST_ELEMENT);
                        newThreeDiagonalMatrix.upperDiagonal.add(buf[0]);
                    }
                    MPI.COMM_WORLD.Recv(buf, 0, 1, MPI.DOUBLE, i, FIRST_FREE_TERM);
                    newFreeTerms.add(buf[0]);
                }
            }

            double newResult[] = new double[newFreeTerms.size()];
            newThreeDiagonalMatrix.size = newFreeTerms.size();
            newThreeDiagonalMatrix.sequentialSolve(newResult, 0,
                    newFreeTerms.stream().mapToDouble(d -> d).toArray());

            reorderParallelResult(newResult);

            int totalIndexDecrease = 0;
            for (int i = 1; i < numOfProcesses; ++i) {
                int iCount = solver.getCount(i) - (i == numOfProcesses - 1 ? 1 : 0);
                MPI.COMM_WORLD.Isend(newResult, 2 * i - 1 - totalIndexDecrease, 2, MPI.DOUBLE, i,
                        UPPER_BOUND_RESULT);

                if (iCount == 1) {
                    ++totalIndexDecrease;
                }

                if (i != numOfProcesses - 1) {
                    MPI.COMM_WORLD.Isend(newResult, 2 * i + 1 - totalIndexDecrease, 2, MPI.DOUBLE, i,
                            LOWER_BOUND_RESULT);
                } else {
                    MPI.COMM_WORLD.Isend(newResult, 2 * i + 1 - totalIndexDecrease, 1, MPI.DOUBLE, i,
                            LOWER_BOUND_RESULT);
                }
            }

            upperBoundResults[0] = newResult[0];
            lowerBoundResults[0] = newResult[1];
            lowerBoundResults[1] = newResult[2];
        }
        else {
            if (count > 1) {
                buf[0] = mainDiagonal.get(0);
                MPI.COMM_WORLD.Isend(buf, 0, 1, MPI.DOUBLE, 0, MAIN_DIAGONAL_FIRST_ELEMENT);
                buf[0] = mainDiagonal.get(size - 1);
                MPI.COMM_WORLD.Isend(buf, 0, 1, MPI.DOUBLE, 0, MAIN_DIAGONAL_LAST_ELEMENT);
                if (rank != numOfProcesses - 1) {
                    MPI.COMM_WORLD.Isend(rightColumn, 0, 1, MPI.DOUBLE, 0, RIGHT_COLUMN_FIRST_ELEMENT);
                    MPI.COMM_WORLD.Isend(rightColumn, rightColumn.length - 1, 1, MPI.DOUBLE, 0,
                            RIGHT_COLUMN_LAST_ELEMENT);
                }
                MPI.COMM_WORLD.Isend(leftColumn, 0, 1, MPI.DOUBLE, 0, LEFT_COLUMN_FIRST_ELEMENT);
                MPI.COMM_WORLD.Isend(leftColumn, leftColumn.length - 1, 1, MPI.DOUBLE, 0,
                        LEFT_COLUMN_LAST_ELEMENT);
                MPI.COMM_WORLD.Isend(freeTerms, 0, 1, MPI.DOUBLE, 0, FIRST_FREE_TERM);
                MPI.COMM_WORLD.Isend(freeTerms, size - 1, 1, MPI.DOUBLE, 0, LAST_FREE_TERM);
            }
            else if (count == 1){
                buf[0] = mainDiagonal.get(0);
                MPI.COMM_WORLD.Isend(buf, 0, 1, MPI.DOUBLE, 0, MAIN_DIAGONAL_FIRST_ELEMENT);
                if (rank != numOfProcesses - 1) {
                    MPI.COMM_WORLD.Isend(rightColumn, 0, 1, MPI.DOUBLE, 0, RIGHT_COLUMN_FIRST_ELEMENT);
                }
                MPI.COMM_WORLD.Isend(leftColumn, 0, 1, MPI.DOUBLE, 0, LEFT_COLUMN_FIRST_ELEMENT);
                MPI.COMM_WORLD.Isend(freeTerms, 0, 1, MPI.DOUBLE, 0, FIRST_FREE_TERM);
            }

            MPI.COMM_WORLD.Recv(upperBoundResults, 0, 2, MPI.DOUBLE, 0,
                    UPPER_BOUND_RESULT);
            if (rank != numOfProcesses - 1) {
                MPI.COMM_WORLD.Recv(lowerBoundResults, 0, 2, MPI.DOUBLE, 0, LOWER_BOUND_RESULT);
            }
            else {
                MPI.COMM_WORLD.Recv(lowerBoundResults, 0, 1, MPI.DOUBLE, 0, LOWER_BOUND_RESULT);
            }
        }

        result[resultOffset] = upperBoundResults[rank == 0 ? 0 : 1];
        result[resultOffset + size - 1] = lowerBoundResults[0];

        if (rank == 0) {
            for (int i = 1; i < size - 1; ++i) {
                result[resultOffset + i] = (freeTerms[i]
                        - rightColumn[i] * lowerBoundResults[1]) / mainDiagonal.get(i - resultOffset);
            }
        }
        else if (rank == numOfProcesses - 1) {
            for (int i = 1; i < size - 1; ++i) {
                result[resultOffset + i] = (freeTerms[i]
                        - leftColumn[i] * upperBoundResults[0]) / mainDiagonal.get(i);
            }
        }
        else {
            for (int i = 1; i < size - 1; ++i) {
                result[resultOffset + i] = (freeTerms[i]
                        - leftColumn[i] * upperBoundResults[0]
                        - rightColumn[i] * lowerBoundResults[1]) / mainDiagonal.get(i);
            }
        }
    }

    private void reorderParallelResult(double[] result) {
        for (int i = 1; i < result.length - 1; i += 2) {
            double temp = result[i];
            result[i] = result[i + 1];
            result[i + 1] = temp;
        }
    }

    @Override
    public void sequentialSolve(double[] result, int resultOffset, double[] freeTerms) {
        if (size != freeTerms.length) {
            throw new RuntimeException("Sizes of matrix and free terms don't match");
        }

        for (int i = 0; i < size - 1; ++i) {
            double alpha = lowerDiagonal.get(i) / mainDiagonal.get(i);
            lowerDiagonal.set(i, 0.0);
            mainDiagonal.set(i + 1, mainDiagonal.get(i + 1) - alpha * upperDiagonal.get(i));
            freeTerms[i + 1] -= alpha * freeTerms[i];
        }

        result[size - 1 + resultOffset] = freeTerms[size - 1] / mainDiagonal.get(size - 1);

        for (int i = size - 2; i >= 0; --i) {
            result[i + resultOffset] = (freeTerms[i] - upperDiagonal.get(i) * result[i + 1 + resultOffset])
                    / mainDiagonal.get(i);
        }
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                stringBuilder.append(get(i, j)).append(' ');
            }
            stringBuilder.append('\n');
        }
        return stringBuilder.toString();
    }
}
