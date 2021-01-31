package diffsolve_test.explicit;

import mpi.MPI;
import numeric_differential_equation.solution.Solution;
import numeric_differential_equation.solver.Solver;

public class TestExplicitParallel {
    public static void main(String[] args) {
        MPI.Init(args);

        Solver solver = MyExplicitEquation.solver;
        long start = System.nanoTime();
        Solution solution = solver.parallel(2).solve();
        long finish = System.nanoTime();
        System.out.println((finish - start) / 1000000);

        MPI.Finalize();
    }
}
