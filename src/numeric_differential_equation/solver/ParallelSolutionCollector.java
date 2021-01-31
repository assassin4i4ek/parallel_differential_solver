package numeric_differential_equation.solver;

import mpi.MPI;
import numeric_differential_equation.solution.ArraySolution;
import numeric_differential_equation.solution.Solution;
import numeric_differential_equation.solution.SolutionArrayImpl;

public interface ParallelSolutionCollector {
    Solution collect(ArraySolution localSolution, Solver solver);

    static ParallelSolutionCollector globalSolutionForAllCollector() {
        return (localSolution, solver) -> {
            SolutionArrayImpl newSolution = new SolutionArrayImpl(solver.mesh, solver.t0, solver.t1, solver.x0, solver.x1);
            for (int i = 0; i < localSolution.getSize(); ++i) {
                int totalCount = 0;
                for (int j = 0; j < solver.counts.length; ++j) {
                    totalCount += solver.counts[j];
                }
                double[] newRow = new double[totalCount];
                double[] row = localSolution.getRow(i);
                MPI.COMM_WORLD.Allgatherv(row, 0, row.length, MPI.DOUBLE, newRow, 0, solver.counts,
                        solver.displacements, MPI.DOUBLE);
                newSolution.setDimRow(newRow, i);
            }

            return newSolution;
        };
    }

    static ParallelSolutionCollector globalSolutionForRankCollector(int rank) {
        return (localSolution, solver) -> {
            SolutionArrayImpl newSolution = null;
            int myRank = MPI.COMM_WORLD.Rank();
            if (myRank == rank) {
                newSolution = new SolutionArrayImpl(solver.mesh, solver.t0, solver.t1, solver.x0, solver.x1);
            }
            for (int i = 0; i < localSolution.getSize(); ++i) {
                int totalCount = 0;
                for (int j = 0; j < solver.counts.length; ++j) {
                    totalCount += solver.counts[j];
                }
                double[] newRow = new double[totalCount];
                double[] row = localSolution.getRow(i);
                MPI.COMM_WORLD.Gatherv(row, 0, row.length, MPI.DOUBLE, newRow, 0, solver.counts,
                        solver.displacements, MPI.DOUBLE, rank);
                if (myRank == rank) {
                    newSolution.setDimRow(newRow, i);
                }
            }

            if (myRank == rank) {
                return newSolution;
            }
            else return null;
        };
    }

    static ParallelSolutionCollector localSolutionCollector() {
        return (localSolution, solver) -> localSolution;
    }

    static ParallelSolutionCollector defaultCollector() {
        return localSolutionCollector();
//        return globalSolutionForRankCollector(0);
//        return globalSolutionForAllCollector();
    }
}
