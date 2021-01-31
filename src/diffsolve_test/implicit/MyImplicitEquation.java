package diffsolve_test.implicit;

import numeric_differential_equation.matrix.Matrix;
import numeric_differential_equation.mesh.Mesh;
import numeric_differential_equation.solution.AccurateSolution;
import numeric_differential_equation.solver.Solver;

public class MyImplicitEquation {
    private static final double t0 = 0;
    private static final double t1 = 1;
    private static final double x0 = 0;
    private static final double x1 = 1;
    private static final double dimStep = 1.0/5;
    private static final double timeStep = 1.0/5;
    static AccurateSolution accurateSolution = (t, x) -> 3 * Math.log(2 * x + 0.5);

    static Solver solver = Solver.builder()
            .implicitNonlinear()
            .withMatrixGenerator(Matrix::threeDiagonalMatrix)
            .forEquation((timeIndex, dimIndex, solution, result) -> {
                double sigma = timeStep / (dimStep * dimStep);
                double leftValue = solution.getValue(timeIndex, dimIndex - 1);
                double rightValue = solution.getValue(timeIndex, dimIndex + 1);
                double middleValue = solution.getValue(timeIndex, dimIndex);
                double prevMiddleValue = solution.getValue(timeIndex - 1, dimIndex);
                double exponentValue = Math.exp(middleValue / 3);
                double freeTerm = -(middleValue - prevMiddleValue
                        - 0.3 * sigma * exponentValue * (rightValue - 2 * middleValue + leftValue
                        + Math.pow(rightValue - leftValue, 2) / 12));
                double leftDerivative = - 0.3 * sigma * exponentValue
                        * (1 - (rightValue - leftValue) / 6);
                double rightDerivative = - 0.3 * sigma * exponentValue
                        * (1 + (rightValue - leftValue) / 6);
                double middleDerivative = 1 - 0.3 * sigma * exponentValue
                        * (-2 + (rightValue - 2 * middleValue + leftValue) / 3
                        + Math.pow(rightValue - leftValue, 2) / 36);
                result.setResult(freeTerm, new double[]{leftDerivative, middleDerivative, rightDerivative});
            }, t0,t1,x0,x1)
            .mesh(Mesh.staticMesh(timeStep, dimStep))
            .initialCondition(x -> 3 * Math.log(2*x + 0.5))
            .boundConditions(t -> 3 * Math.log(0.5), t -> 3 * Math.log(2.5));
}
