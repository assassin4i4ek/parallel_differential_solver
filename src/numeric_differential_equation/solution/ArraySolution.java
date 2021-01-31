package numeric_differential_equation.solution;

public interface ArraySolution extends Solution {
    double[] getRow(int timeIndex);
}
