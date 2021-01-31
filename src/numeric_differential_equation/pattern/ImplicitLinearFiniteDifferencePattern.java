package numeric_differential_equation.pattern;

import numeric_differential_equation.pattern.pattern_result.ImplicitLinearPatternResult;
import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.solution.Solution;

public interface ImplicitLinearFiniteDifferencePattern extends FiniteDifferencePattern {
    void applyPattern(int timeIndex, int dimIndex, Solution solution, ImplicitLinearPatternResult result);

    @Override
    default void applyPattern(int timeIndex, int dimIndex, Solution solution, PatternResult result) {
        applyPattern(timeIndex, dimIndex, solution, (ImplicitLinearPatternResult) result);
    }
}
