package diffsolve_test.implicit;

import mpi.MPI;
import numeric_differential_equation.solution.Solution;
import numeric_differential_equation.solver.Solver;

public class TestImplicitParallel {
    public static void main(String[] args) {
        MPI.Init(args);

        Solver solver = MyImplicitEquation.solver;
        long start = System.nanoTime();
        Solution solution = solver
                .parallel(4)
                .solve();
//        if (MPI.COMM_WORLD.Rank() == 0) {
//            System.out.println(solution);
//        }
        long finish = System.nanoTime();
        System.out.println((finish - start) / 1000000);

        MPI.Finalize();
    }
}
