package numeric_differential_equation.pattern;

import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.solution.Solution;

public interface FiniteDifferencePattern {
    void applyPattern(int timeIndex, int dimIndex, Solution solution, PatternResult result);
}
