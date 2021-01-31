package numeric_differential_equation.solver;

import numeric_differential_equation.matrix.Matrix;
import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.solution.SolutionArrayImpl;

import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinWorkerThread;
import java.util.function.Supplier;
import java.util.stream.IntStream;

public class ImplicitLinearSolver extends Solver {
    Supplier<Matrix> matrixGenerator;

    @Override
    protected void sequentialSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps, SolutionArrayImpl solution) {
        double solutionValues[] = new double[dimSteps.size()];
        double freeTerms[] = new double[dimSteps.size() - 2];
        int dimIndex;
        double leftBoundValue = leftBoundTimeCondition.apply(t);
        double rightBoundValue = rightBoundTimeCondition.apply(t);
        solutionValues[0] = leftBoundValue;
        solutionValues[solutionValues.length - 1] = rightBoundValue;
        Matrix matrix = matrixGenerator.get();

        pattern.applyPattern(timeIndex, 1, solution, patternResults[0]);
        matrix.setRow(leftShift(patternResults[0].getCoefficients()), 0, 0);
        freeTerms[0] = patternResults[0].getFreeTerm() - leftBoundValue * patternResults[0].getCoefficients()[0];

        for (dimIndex = 2; dimIndex < dimSteps.size() - 2; ++dimIndex) {
            pattern.applyPattern(timeIndex, dimIndex, solution, patternResults[0]);
            matrix.setRow(patternResults[0].getCoefficients(), dimIndex - 1, dimIndex - 2);
            freeTerms[dimIndex - 1] = patternResults[0].getFreeTerm();
        }

        pattern.applyPattern(timeIndex, dimIndex, solution, patternResults[0]);
        matrix.setRow(rightShift(patternResults[0].getCoefficients()), dimIndex - 1, dimIndex - 2);
        freeTerms[freeTerms.length - 1] = patternResults[0].getFreeTerm() - rightBoundValue * patternResults[0]
                .getCoefficients()[patternResults[0].getCoefficients().length - 1];

        matrix.sequentialSolve(solutionValues, 1, freeTerms);
        solution.setDimRow(solutionValues, timeIndex);
    }

    @Override
    protected void parallelSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps, SolutionArrayImpl solution) {
        double solutionValues[] = new double[getCount(rank)];
        Matrix matrix = matrixGenerator.get();

        int fromIndex = myFromIndex + (rank == 0 ? 1 : 0);
        int toIndex = myToIndex - (rank == numOfProcesses - 1 ? 1 : 0);
        int dimIndexFrom = fromIndex + (rank == 0 ? 1 : 0);
        int dimIndexTo = toIndex - (rank == numOfProcesses - 1 ? 1 : 0);
        double freeTerms[] = new double[toIndex - fromIndex];
        matrix.initSize(freeTerms.length);

        if (rank == 0) {
            double leftBoundValue = leftBoundTimeCondition.apply(t);

            solutionValues[0] = leftBoundValue;
            pattern.applyPattern(timeIndex, fromIndex, solution, patternResults[0]);
            matrix.setRow(leftShift(patternResults[0].getCoefficients()), 0, 0);
            freeTerms[0] = patternResults[0].getFreeTerm() - leftBoundValue * patternResults[0].getCoefficients()[0];
        }

        try {
            forkJoinPool.submit(() -> IntStream.range(dimIndexFrom, dimIndexTo)
                    .parallel()
                    .forEach(dimIndex -> {
                        PatternResult patternResult = patternResults[((ForkJoinWorkerThread) Thread.currentThread())
                                .getPoolIndex()];
                        pattern.applyPattern(timeIndex, dimIndex, solution, patternResult);
                        matrix.setRow(patternResult.getCoefficients(), dimIndex - fromIndex,
                                dimIndex - fromIndex - 1);
                        freeTerms[dimIndex - fromIndex] = patternResult.getFreeTerm();
                    })).get();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        if (rank == numOfProcesses - 1) {
            double rightBoundValue = rightBoundTimeCondition.apply(t);

            solutionValues[solutionValues.length - 1] = rightBoundValue;
            pattern.applyPattern(timeIndex, toIndex - 1, solution, patternResults[0]);
            matrix.setRow(rightShift(patternResults[0].getCoefficients()), toIndex - fromIndex - 1,
                    toIndex - fromIndex - 2);
            freeTerms[freeTerms.length - 1] = patternResults[0].getFreeTerm() - rightBoundValue
                    * patternResults[0].getCoefficients()[patternResults[0].getCoefficients().length - 1];
        }

        matrix.parallelSolve(solutionValues, rank == 0 ? 1 : 0, freeTerms, this);
        solution.setDimRow(solutionValues, timeIndex);
    }
}
