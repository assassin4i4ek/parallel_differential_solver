package numeric_differential_equation.solver;

import mpi.MPI;
import numeric_differential_equation.matrix.Matrix;
import numeric_differential_equation.mpi_features.MpiTags;
import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.solution.SolutionArrayImpl;

import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ForkJoinWorkerThread;
import java.util.function.Supplier;
import java.util.stream.IntStream;

public class ImplicitNonlinearSolver extends Solver {
    Supplier<Matrix> matrixGenerator;

    @Override
    protected void sequentialSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps, SolutionArrayImpl solution) {
        double solutionValues[] = Arrays.copyOf(solution.getRow(timeIndex - 1), dimSteps.size());//new double[dimSteps.size()];
        int dimIndex;
        double leftBoundValue = leftBoundTimeCondition.apply(t);
        double rightBoundValue = rightBoundTimeCondition.apply(t);
        solutionValues[0] = leftBoundValue;
        solutionValues[solutionValues.length - 1] = rightBoundValue;
        solution.setDimRow(solutionValues, timeIndex);

        double freeTerms[] = new double[dimSteps.size()];
        double deltaSolutionValues[] = new double[dimSteps.size()];

        do {
            Matrix matrix = matrixGenerator.get();
            matrix.initSize(freeTerms.length);
            matrix.set(0,0,1);
            for (dimIndex = 1; dimIndex < dimSteps.size() - 1; ++dimIndex) {
                pattern.applyPattern(timeIndex, dimIndex, solution, patternResults[0]);
                matrix.setRow(patternResults[0].getCoefficients(), dimIndex, dimIndex - 1);
                freeTerms[dimIndex] = patternResults[0].getFreeTerm();
            }

            matrix.set(freeTerms.length - 1, freeTerms.length - 1, 1);

            matrix.sequentialSolve(deltaSolutionValues, 0, freeTerms);

            for (int i = 0; i < solutionValues.length; ++i) {
                solutionValues[i] += deltaSolutionValues[i];
            }
        }
        while (!sequentialArrayApproximatelyZero(deltaSolutionValues));
    }

    @Override
    protected void parallelSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps, SolutionArrayImpl solution) {
        int fromIndex = myFromIndex + (rank == 0 ? 1 : 0);
        int toIndex = myToIndex - (rank == numOfProcesses - 1 ? 1 : 0);
        double freeTerms[] = new double[myToIndex - myFromIndex];
        double deltaSolutionValues[] = new double[myToIndex - myFromIndex];
        double solutionValues[] = Arrays.copyOf(solution.getRow(timeIndex - 1), myToIndex - myFromIndex);
        if (rank == 0) {
            solutionValues[0] = leftBoundTimeCondition.apply(t);
        }

        if (rank == numOfProcesses - 1) {
            solutionValues[solutionValues.length - 1] = rightBoundTimeCondition.apply(t);
        }

        solution.setDimRow(solutionValues, timeIndex);

        do {
            Matrix matrix = matrixGenerator.get();
            matrix.initSize(freeTerms.length);
            if (rank == 0) {
                matrix.set(0,0,1);
            }

            forkJoinPool.submit(() -> IntStream.range(fromIndex, toIndex)
                    .parallel()
                    .forEach(dimIndex -> {
                        PatternResult patternResult = patternResults[((ForkJoinWorkerThread) Thread.currentThread())
                                .getPoolIndex()];
                        pattern.applyPattern(timeIndex, dimIndex, solution, patternResult);
                        matrix.setRow(patternResult.getCoefficients(), dimIndex - myFromIndex,
                                dimIndex - myFromIndex - 1);
                        freeTerms[dimIndex - myFromIndex] = patternResult.getFreeTerm();
                    })).join();

            if (rank == numOfProcesses - 1) {
                matrix.set(freeTerms.length - 1, freeTerms.length - 1, 1);
            }

            matrix.parallelSolve(deltaSolutionValues, 0, freeTerms, this);
            for (int i = 0; i < deltaSolutionValues.length; ++i) {
                solutionValues[i] += deltaSolutionValues[i];
            }
            System.out.println();

            MPI.COMM_WORLD.Barrier();
        }
        while (!parallelArrayApproximatelyZero(deltaSolutionValues));
    }

    private boolean parallelArrayApproximatelyZero(double[] array) {
        boolean[] result = {sequentialArrayApproximatelyZero(array)};
        if (rank == 0) {
            for (int i = 1; i < numOfProcesses && result[0]; ++i) {
                MPI.COMM_WORLD.Recv(result, 0,1,MPI.BOOLEAN, i, MpiTags.ARRAY_APPROXIMATELY_ZERO);
            }

            for (int i = 1; i < numOfProcesses; ++i) {
                MPI.COMM_WORLD.Isend(result, 0,1,MPI.BOOLEAN, i, MpiTags.ARRAY_APPROXIMATELY_ZERO);
            }
        }
        else {
            MPI.COMM_WORLD.Isend(result, 0,1,MPI.BOOLEAN,0,MpiTags.ARRAY_APPROXIMATELY_ZERO);
            MPI.COMM_WORLD.Recv(result,0,1, MPI.BOOLEAN, 0, MpiTags.ARRAY_APPROXIMATELY_ZERO);
        }
        return result[0];
    }

    private boolean sequentialArrayApproximatelyZero(double[] array) {
        for (double item : array) {
            if (Math.abs(item) > 0.00000001) {
                return false;
            }
        }

        return true;
    }
}
