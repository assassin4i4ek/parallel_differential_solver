package numeric_differential_equation.pattern.pattern_result;

public class PatternResultImpl implements PatternResult {
    private double result;
    private double[] coefficients;

    @Override
    public void setResult(double value) {
        result = value;
    }

    @Override
    public double getResult() {
        return result;
    }

    @Override
    public void setResult(double value, double [] coefficients) {
        result = value;
        this.coefficients = coefficients;
    }

    @Override
    public double[] getCoefficients() {
        return coefficients;
    }

    @Override
    public double getFreeTerm() {
        return result;
    }
}
