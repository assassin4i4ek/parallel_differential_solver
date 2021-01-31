package diffsolve_test.explicit;

import numeric_differential_equation.solution.AccurateSolution;
import numeric_differential_equation.solution.Solution;
import numeric_differential_equation.solver.Solver;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class TestExplicitSequential {
    public static void main(String[] args) {
        Solver solver = MyExplicitEquation.solver;
        AccurateSolution accurateSolution = MyExplicitEquation.accurateSolution;
        long start = System.nanoTime();
        Solution solution = solver.sequential().solve();
        long finish = System.nanoTime();
        System.out.println((finish - start) / 1000000);
//        System.out.println(solution);
//        solution.printError(accurateSolution);

        try {
            PrintWriter writer = new PrintWriter(new FileWriter("C:\\Users\\Admin\\Desktop\\ParComp\\Course_new\\src\\1.txt"));
            PrintWriter writer2 = new PrintWriter(new FileWriter("C:\\Users\\Admin\\Desktop\\ParComp\\Course_new\\src\\2.txt"));
            solution.printSolution(writer);
            solution.printAccurateSolution(accurateSolution, writer2);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
//        System.out.println("Average absolute error: " + String.format("%.8f", solution.meanAbsError(accurateSolution)));
//        System.out.println("Maximum absolute error: " + String.format("%.8f", solution.maxAbsError(accurateSolution)));
//        System.out.println("Average relative error: " + String.format("%.8f", solution.meanRelativeError(accurateSolution)));
//        System.out.println("Maximum relative error: " + String.format("%.8f", solution.maxRelativeError(accurateSolution)));
    }
}
