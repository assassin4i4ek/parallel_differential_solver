package diffsolve_test.explicit;

import numeric_differential_equation.mesh.Mesh;
import numeric_differential_equation.solution.AccurateSolution;
import numeric_differential_equation.solver.Solver;

public class MyExplicitEquation {
    private static final double t0 = 0;
    private static final double t1 = 1;
    private static final double x0 = 0;
    private static final double x1 = 1;
    private static final double dimStep = 1.0/5;
    private static final double timeStep = 1.0/50;
    static AccurateSolution accurateSolution = (t, x) -> 3 * Math.log(2 * x + 0.5);

    static Solver solver = Solver.builder()
            .explicit()
            .forEquation((timeIndex, dimIndex, solution, result) -> {
                double sigma = timeStep / (dimStep * dimStep);
                double leftValue = solution.getValue(timeIndex - 1, dimIndex - 1);
                double rightValue = solution.getValue(timeIndex - 1, dimIndex + 1);
                double middleValue = solution.getValue(timeIndex - 1, dimIndex);
                result.setResult(middleValue
                        + 0.3 * sigma * Math.exp(middleValue / 3)
                        * (rightValue - 2 * middleValue + leftValue + Math.pow(rightValue - leftValue, 2)/12));
            }, t0,t1,x0,x1)
            .mesh(Mesh.staticMesh(timeStep, dimStep))
            .initialCondition(x -> 3 * Math.log(2*x + 0.5))
            .boundConditions(t -> 3 * Math.log(0.5), t -> 3 * Math.log(2.5));
}
