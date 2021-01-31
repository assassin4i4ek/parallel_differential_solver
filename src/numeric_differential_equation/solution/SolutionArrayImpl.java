package numeric_differential_equation.solution;

import numeric_differential_equation.DoubleArray;
import numeric_differential_equation.mesh.Mesh;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class SolutionArrayImpl implements ArraySolution {
    protected List<DoubleArray> matrix = new ArrayList<>();
    private Mesh mesh;
    private double t0;
    private double t1;
    private double x0;
    private double x1;

    public SolutionArrayImpl(Mesh mesh, double t0, double t1, double x0, double x1) {
        this.mesh = mesh;
        this.t0 = t0;
        this.t1 = t1;
        this.x0 = x0;
        this.x1 = x1;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof Solution) {
            for (int i = 0, matrixSize = matrix.size(); i < matrixSize; i++) {
                DoubleArray doubleArray = matrix.get(i);
                double[] array = doubleArray.array;
                for (int j = 0; j < array.length; j++) {
                    double d = array[j];
                    if (Math.abs(d - ((Solution) obj).getValue(i,j)) > 0.000001) {
                        return false;
                    }
                }
            }
            return true;
        }
        return false;
    }

    @Override
    public double getValue(int timeIndex, int dimIndex) {
        return matrix.get(timeIndex).array[dimIndex];
    }

    @Override
    public int getSize() {
        return matrix.size();
    }

    @Override
    public void printError(AccurateSolution accurateSolution) {
        double t = t0;
        for (int timeIndex = 0; timeIndex < matrix.size(); ++timeIndex, t = mesh.nextTimeStep(t)) {
            double x = x0;
            for (int dimIndex = 0; dimIndex < matrix.get(timeIndex).array.length; ++dimIndex, x = mesh.nextDimStep(x, t)) {
                double accurateValue = accurateSolution.apply(t, x);
                double solutionValue = getValue(timeIndex, dimIndex);
                double localError = Math.abs(solutionValue - accurateValue);
                System.out.print(String.format("%.4f", localError) + " ");
            }
            System.out.println();
        }
    }

    @Override
    public double meanAbsError(AccurateSolution accurateSolution) {
        double absError = 0;
        long counter = 0;
        double t = t0;
        for (int timeIndex = 0; timeIndex < matrix.size(); ++timeIndex, t = mesh.nextTimeStep(t)) {
            double x = x0;
            for (int dimIndex = 0; dimIndex < matrix.get(timeIndex).array.length; ++dimIndex, x = mesh.nextDimStep(x, t)) {
                double accurateValue = accurateSolution.apply(t, x);
                double solutionValue = getValue(timeIndex, dimIndex);
                double localError = Math.abs(solutionValue - accurateValue);
                absError = (absError / (counter + 1)) * counter + localError / (counter + 1);
                ++counter;
            }
        }

        return absError;
    }

    @Override
    public double maxAbsError(AccurateSolution accurateSolution) {
        double maxError = Double.MIN_VALUE;
        double t = t0;
        for (int timeIndex = 0; timeIndex < matrix.size(); ++timeIndex, t = mesh.nextTimeStep(t)) {
            double x = x0;
            for (int dimIndex = 0; dimIndex < matrix.get(timeIndex).array.length; ++dimIndex, x = mesh.nextDimStep(x, t)) {
                double accurateValue = accurateSolution.apply(t, x);
                double solutionValue = getValue(timeIndex, dimIndex);
                double localError = Math.abs(solutionValue - accurateValue);
                if (localError > maxError) {
                    maxError = localError;
                }
            }
        }
        return maxError;
    }

    @Override
    public double meanRelativeError(AccurateSolution accurateSolution) {
        double relError = 0;
        long counter = 0;
        double t = t0;
        for (int timeIndex = 0; timeIndex < matrix.size(); ++timeIndex, t = mesh.nextTimeStep(t)) {
            double x = x0;
            for (int dimIndex = 0; dimIndex < matrix.get(timeIndex).array.length; ++dimIndex, x = mesh.nextDimStep(x, t)) {
                double accurateValue = accurateSolution.apply(t, x);
                double solutionValue = getValue(timeIndex, dimIndex);
                if (accurateValue != 0.0 && solutionValue != 0.0) {
                    double localError = Math.abs((solutionValue - accurateValue) / accurateValue);
                    relError += localError;
                    ++counter;
                }
            }
        }
        relError /= counter;
        return relError;
    }

    @Override
    public double maxRelativeError(AccurateSolution accurateSolution) {
        double maxError = Double.MIN_VALUE;
        double t = t0;
        for (int timeIndex = 0; timeIndex < matrix.size(); ++timeIndex, t = mesh.nextTimeStep(t)) {
            double x = x0;
            for (int dimIndex = 0; dimIndex < matrix.get(timeIndex).array.length; ++dimIndex, x = mesh.nextDimStep(x, t)) {
                double accurateValue = accurateSolution.apply(t, x);
                double solutionValue = getValue(timeIndex, dimIndex);
                if (accurateValue != 0.0 && solutionValue != 0.0) {
                    double localError = Math.abs(solutionValue - accurateValue) / accurateValue;
                    if (localError > maxError) {
                        maxError = localError;
                    }
                }
            }
        }
        return maxError;
    }

    @Override
    public double[] getRow(int timeIndex) {
        return matrix.get(timeIndex).array;
    }

    public void setDimRow(double[] dimRow, int timeIndex) {
        for (int k = matrix.size(); k < timeIndex; ++k) {
            matrix.add(null);
        }
        matrix.add(timeIndex, new DoubleArray(dimRow));
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        for (DoubleArray list : matrix) {
            for (double value : list.array) {
                builder.append(String.format("%.4f", value)).append("\t");
            }
            builder.append('\n');
        }
        return builder.toString();
    }

    @Override
    public void printSolution(PrintWriter writer) {
        writer.write("[");
        for (DoubleArray list : matrix) {
            double[] array = list.array;
            for (double value : array) {
                writer.print(value);
                writer.print(' ');
            }
            writer.print(";");
        }
        writer.write("]");
        writer.flush();
    }

    @Override
    public void printAccurateSolution(AccurateSolution accurateSolution, PrintWriter writer) {
        writer.write("[");

        double t = t0;
        for (int timeIndex = 0; timeIndex < matrix.size(); ++timeIndex, t = mesh.nextTimeStep(t)) {
            double x = x0;
            for (int dimIndex = 0; dimIndex < matrix.get(timeIndex).array.length; ++dimIndex, x = mesh.nextDimStep(x, t)) {
                double accurateValue = accurateSolution.apply(t, x);
                writer.print(accurateValue);
                writer.print(' ');
            }
            writer.print(";");
        }

        writer.write("]");
        writer.flush();
    }
}
