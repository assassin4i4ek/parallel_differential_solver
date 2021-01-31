package numeric_differential_equation.solver;

import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.solution.SolutionArrayImpl;

import java.util.List;
import java.util.concurrent.ForkJoinWorkerThread;
import java.util.stream.IntStream;

public class ExplicitSolver extends Solver {
    @Override
    protected void sequentialSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps, SolutionArrayImpl solution) {
        double solutionValues[] = new double[dimSteps.size()];

        solutionValues[0] = leftBoundTimeCondition.apply(t);
        for (int dimIndex = 1; dimIndex < solutionValues.length - 1; ++dimIndex) {
            pattern.applyPattern(timeIndex, dimIndex, solution, patternResults[0]);
            solutionValues[dimIndex] = patternResults[0].getResult();
        }
        solutionValues[solutionValues.length - 1] = rightBoundTimeCondition.apply(t);

        solution.setDimRow(solutionValues, timeIndex);
    }

    @Override
    protected void parallelSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps, SolutionArrayImpl solution) {
        double solutionValues[] = new double[getCount(rank)];
        if (rank == 0) {
            solutionValues[0] = leftBoundTimeCondition.apply(t);
        }
        else if (rank == numOfProcesses - 1) {
            solutionValues[solutionValues.length - 1] = rightBoundTimeCondition.apply(t);
        }

        int fromIndex = myFromIndex + (rank == 0 ? 1 : 0);
        int toIndex = myToIndex - (rank == numOfProcesses - 1 ? 1 : 0);

        forkJoinPool.submit(() -> IntStream.range(fromIndex, toIndex)
                    .parallel()
                    .forEach(dimIndex -> {
                        PatternResult patternResult = patternResults[((ForkJoinWorkerThread) Thread.currentThread())
                                .getPoolIndex()];
                        pattern.applyPattern(timeIndex, dimIndex, solution, patternResult);
                        solutionValues[dimIndex - myFromIndex] = patternResult.getResult();
                    })).join();

        solution.setDimRow(solutionValues, timeIndex);
    }
}
