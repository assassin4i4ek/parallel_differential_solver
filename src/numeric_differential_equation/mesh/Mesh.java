package numeric_differential_equation.mesh;

public interface Mesh {
    double nextTimeStep(double prevTimeStep);
    double nextDimStep(double prevDimStep, double curTimeStep);

    static StaticMesh staticMesh(double timeStep, double dimStep) {
        return new StaticMesh() {
            @Override
            public double nextTimeStep(double prevTimeStep) {
                return prevTimeStep + timeStep;
            }

            @Override
            public double nextDimStep(double prevDimStep, double curTimeStep) {
                return prevDimStep + dimStep;
            }
        };
    }
}
