package numeric_differential_equation.pattern.pattern_result;

public interface PatternResult extends ImplicitLinearPatternResult, ExplicitPatternResult,
        ImplicitNonlinearPatternResult {
    double getResult();
    double[] getCoefficients();
    double getFreeTerm();
}
