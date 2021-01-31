package numeric_differential_equation.pattern;

import numeric_differential_equation.pattern.pattern_result.ImplicitNonlinearPatternResult;
import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.solution.Solution;

public interface ImplicitNonlinearFiniteDifferencePattern extends FiniteDifferencePattern {
    void applyPattern(int timeIndex, int dimIndex, Solution solution, ImplicitNonlinearPatternResult result);

    @Override
    default void applyPattern(int timeIndex, int dimIndex, Solution solution, PatternResult result) {
        applyPattern(timeIndex, dimIndex, solution, (ImplicitNonlinearPatternResult) result);
    }
}
