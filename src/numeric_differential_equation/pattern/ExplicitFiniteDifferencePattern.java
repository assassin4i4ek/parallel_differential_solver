package numeric_differential_equation.pattern;

import numeric_differential_equation.pattern.pattern_result.ExplicitPatternResult;
import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.solution.Solution;

public interface ExplicitFiniteDifferencePattern extends FiniteDifferencePattern {
    void applyPattern(int timeIndex, int dimIndex, Solution solution, ExplicitPatternResult result);

    @Override
    default void applyPattern(int timeIndex, int dimIndex, Solution solution, PatternResult result) {
        applyPattern(timeIndex, dimIndex, solution, (ExplicitPatternResult) result);
    }
}
