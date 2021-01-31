package numeric_differential_equation.solution;

import mpi.MPI;
import numeric_differential_equation.mesh.Mesh;
import numeric_differential_equation.mpi_features.MpiTags;
import numeric_differential_equation.solver.Solver;

public class ParallelSolutionArrayImpl extends SolutionArrayImpl {
    private Solver solver;
    private int fromIndex;
    private int toIndex;

    public ParallelSolutionArrayImpl(Solver solver, Mesh mesh, double t0, double t1, double x0, double x1) {
        super(mesh, t0, t1, x0, x1);
        this.solver = solver;
        fromIndex = solver.getFromIndex();
        toIndex = solver.getToIndex();
    }

    @Override
    public double getValue(int timeIndex, int dimIndex) {
        if (fromIndex <= dimIndex && dimIndex < toIndex) {
            while (timeIndex >= matrix.size()) {
                synchronized (this) {
                    while (timeIndex >= matrix.size()) {
                        try {
                            wait();
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    }
                }
            }
            return super.getValue(timeIndex, dimIndex - fromIndex);
        }
        else {
            int[] outBuf = {timeIndex, dimIndex};
            MPI.COMM_WORLD.Isend(outBuf, 0, 2, MPI.INT, solver.getRankForDimIndex(dimIndex),
                    MpiTags.SOLUTION_VALUE_REQUEST);
            double[] inBuf = new double[1];
            MPI.COMM_WORLD.Recv(inBuf, 0, 1, MPI.DOUBLE, solver.getRankForDimIndex(dimIndex),
                    MpiTags.SOLUTION_VALUE_REPLY);
            return inBuf[0];
        }
    }

    @Override
    public void setDimRow(double[] dimRow, int timeIndex) {
        super.setDimRow(dimRow, timeIndex);
        synchronized (this) {
            notifyAll();
        }
    }

    @Override
    public double meanAbsError(AccurateSolution accurateSolution) {
        return super.meanAbsError(accurateSolution);
    }
}
