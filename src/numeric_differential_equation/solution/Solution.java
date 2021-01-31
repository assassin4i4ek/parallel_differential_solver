package numeric_differential_equation.solution;

import java.io.PrintWriter;

public interface Solution {
    double getValue(int timeIndex, int dimIndex);
    int getSize();
    double meanAbsError(AccurateSolution accurateSolution);
    double maxAbsError(AccurateSolution accurateSolution);
    double meanRelativeError(AccurateSolution accurateSolution);
    double maxRelativeError(AccurateSolution accurateSolution);
    void printError(AccurateSolution accurateSolution);
    void printSolution(PrintWriter writer);
    void printAccurateSolution(AccurateSolution accurateSolution, PrintWriter writer);
}
