package numeric_differential_equation.solver;

import mpi.MPI;
import mpi.Status;
import numeric_differential_equation.functions.OneDimensionFunction;
import numeric_differential_equation.matrix.Matrix;
import numeric_differential_equation.mesh.Mesh;
import numeric_differential_equation.mesh.StaticMesh;
import numeric_differential_equation.mpi_features.MpiTags;
import numeric_differential_equation.pattern.ExplicitFiniteDifferencePattern;
import numeric_differential_equation.pattern.FiniteDifferencePattern;
import numeric_differential_equation.pattern.ImplicitLinearFiniteDifferencePattern;
import numeric_differential_equation.pattern.ImplicitNonlinearFiniteDifferencePattern;
import numeric_differential_equation.pattern.pattern_result.PatternResult;
import numeric_differential_equation.pattern.pattern_result.PatternResultImpl;
import numeric_differential_equation.solution.ParallelSolutionArrayImpl;
import numeric_differential_equation.solution.Solution;
import numeric_differential_equation.solution.SolutionArrayImpl;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.function.Supplier;
import java.util.stream.IntStream;

public abstract class Solver {
    double t0;
    double t1;
    double x0;
    double x1;
    FiniteDifferencePattern pattern;
    private boolean parallelismParam;
    Mesh mesh;
    private OneDimensionFunction initialDimensionCondition;
    OneDimensionFunction leftBoundTimeCondition;
    OneDimensionFunction rightBoundTimeCondition;
    PatternResult patternResults[];
    //parallelism
    private int forkJoinPoolSize;
    ForkJoinPool forkJoinPool;
    private ParallelSolutionCollector parallelSolutionCollector = null;
    int numOfProcesses;
    int rank = 0;
    int[] counts = new int[0];
    int[] displacements = new int[0];
    int myFromIndex = 0;
    int myToIndex = 0;
    private Thread supportThread;

    public int getFromIndex() {
        return myFromIndex;
    }

    public int getToIndex() {
        return myToIndex;
    }

    public int getNumOfProcesses() {
        return numOfProcesses;
    }

    public int getRank() {
        return rank;
    }

    public int getCount(int rank) {
        return counts[rank];
    }

    public int getRankForDimIndex(int dimIndex) {
        int iRank;
        for (iRank = 0; iRank < numOfProcesses - 1 && dimIndex >= displacements[iRank + 1]; ++iRank);
        return iRank;
    }

    public interface ImplicitnessBuilder {
        ImplicitLinearMatrixBuilder implicitLinear();
        ImplicitNonlinearMatrixBuilder implicitNonlinear();
        ExplicitEquationBuilder explicit();
    }

    public interface ImplicitNonlinearMatrixBuilder {
        ImplicitNonlinearEquationBuilder withMatrixGenerator(Supplier<Matrix> matrixGenerator);
    }

    public interface ImplicitLinearMatrixBuilder {
        ImplicitLinearEquationBuilder withMatrixGenerator(Supplier<Matrix> matrixGenerator);
    }

    public interface ImplicitLinearEquationBuilder {
        MeshBuilder forEquation(ImplicitLinearFiniteDifferencePattern pattern,
                                double t0, double t1, double x0, double x1);
    }

    public interface ImplicitNonlinearEquationBuilder {
        MeshBuilder forEquation(ImplicitNonlinearFiniteDifferencePattern pattern,
                                double t0, double t1, double x0, double x1);
    }

    public interface ExplicitEquationBuilder {
        MeshBuilder forEquation(ExplicitFiniteDifferencePattern pattern, double t0, double t1, double x0, double x1);
    }

    public interface MeshBuilder {
        InitialConditionBuilder mesh(Mesh mesh);
    }

    public interface InitialConditionBuilder {
        BoundConditionsBuilder initialCondition(OneDimensionFunction initialDimensionCondition);
    }

    public interface BoundConditionsBuilder {
        Solver boundConditions(OneDimensionFunction leftBoundCondition,
                               OneDimensionFunction rightBoundCondition);
    }

    public static ImplicitnessBuilder builder() {
        Solver instance[] = new Solver[1];

        BoundConditionsBuilder boundConditionsBuilder = (leftBoundCondition, rightBoundCondition) -> {
            instance[0].leftBoundTimeCondition = leftBoundCondition;
            instance[0].rightBoundTimeCondition = rightBoundCondition;
            return instance[0];
        };
        InitialConditionBuilder initialConditionBuilder = initialDimensionCondition -> {
            instance[0].initialDimensionCondition = initialDimensionCondition;
            return boundConditionsBuilder;
        };
        MeshBuilder meshBuilder = mesh -> {
            instance[0].mesh = mesh;
            return initialConditionBuilder;
        };
        ImplicitLinearEquationBuilder implicitLinearEquationBuilder = (pattern, t0, t1, x0, x1) -> {
            instance[0].pattern = pattern;
            instance[0].t0 = t0;
            instance[0].t1 = t1;
            instance[0].x0 = x0;
            instance[0].x1 = x1;
            return meshBuilder;
        };
        ImplicitNonlinearEquationBuilder implicitNonlinearEquationBuilder = (pattern, t0, t1, x0, x1) -> {
            instance[0].pattern = pattern;
            instance[0].t0 = t0;
            instance[0].t1 = t1;
            instance[0].x0 = x0;
            instance[0].x1 = x1;
            return meshBuilder;
        };
        ExplicitEquationBuilder explicitEquationBuilder = (pattern, t0, t1, x0, x1) -> {
            instance[0].pattern = pattern;
            instance[0].t0 = t0;
            instance[0].t1 = t1;
            instance[0].x0 = x0;
            instance[0].x1 = x1;
            return meshBuilder;
        };

        ImplicitLinearMatrixBuilder linearMatrixBuilder = matrixGenerator -> {
            ((ImplicitLinearSolver) instance[0]).matrixGenerator = matrixGenerator;
            return implicitLinearEquationBuilder;
        };

        ImplicitNonlinearMatrixBuilder nonlinearMatrixBuilder = matrixGenerator -> {
            ((ImplicitNonlinearSolver) instance[0]).matrixGenerator = matrixGenerator;
            return implicitNonlinearEquationBuilder;
        };

        return new ImplicitnessBuilder() {
            @Override
            public ImplicitLinearMatrixBuilder implicitLinear() {
                instance[0] = new ImplicitLinearSolver();
                return linearMatrixBuilder;
            }

            @Override
            public ImplicitNonlinearMatrixBuilder implicitNonlinear() {
                instance[0] = new ImplicitNonlinearSolver();
                return nonlinearMatrixBuilder;
            }

            @Override
            public ExplicitEquationBuilder explicit() {
                instance[0] = new ExplicitSolver();
                return explicitEquationBuilder;
            }
        };
    }

    public Solver sequential() {
        parallelismParam = false;
        return this;
    }

    public Solver parallel() {
        forkJoinPoolSize = Runtime.getRuntime().availableProcessors();
        parallelismParam = true;
        return this;
    }

    public Solver parallel(ParallelSolutionCollector collector) {
        parallelSolutionCollector = collector;
        forkJoinPoolSize = Runtime.getRuntime().availableProcessors();
        return this;
    }

    public Solver parallel(int forkJoinPoolSize) {
        this.forkJoinPoolSize = forkJoinPoolSize;
        parallelismParam = true;
        return this;
    }

    public Solver parallel(int forkJoinPoolSize, ParallelSolutionCollector collector) {
        parallelSolutionCollector = collector;
        this.forkJoinPoolSize = forkJoinPoolSize;
        parallelismParam = true;
        return this;
    }

    public Solution solve() {
        SolutionArrayImpl solution;
        double t = mesh.nextTimeStep(t0);
        int timeIndex = 1;
        List<Double> dimSteps = getDimSteps(t0);
        //MPI params
        if (parallelismParam) {
            forkJoinPool = new ForkJoinPool(forkJoinPoolSize);
            patternResults = new PatternResult[forkJoinPool.getParallelism()];
            for (int i = 0; i < patternResults.length; ++i) {
                patternResults[i] = new PatternResultImpl();
            }
            //init MPI features
            numOfProcesses = MPI.COMM_WORLD.Size();
            rank = MPI.COMM_WORLD.Rank();
            counts = new int[numOfProcesses];
            displacements = new int[numOfProcesses];
            initMpi(dimSteps.size());
            myFromIndex = displacements[rank];
            myToIndex = myFromIndex + counts[rank];
            //solution
            solution = new ParallelSolutionArrayImpl(this, mesh, t0, t1, x0, x1);
            if (parallelSolutionCollector == null) {
                parallelSolutionCollector = ParallelSolutionCollector.defaultCollector();
            }

            //support thread
            if (supportThread != null && supportThread.isAlive()) {
                supportThread.interrupt();
            }
            supportThread = new Thread(() -> {
                ThreadPoolExecutor supportPool = (ThreadPoolExecutor) Executors.newCachedThreadPool((runnable) -> {
                    Thread thread = Executors.defaultThreadFactory().newThread(runnable);
                    thread.setDaemon(true);
                    return thread;
                });

                int inBuf[] = new int[2];
                while (true) {
                    Status status = MPI.COMM_WORLD.Recv(inBuf, 0, 2, MPI.INT, MPI.ANY_SOURCE,
                            MpiTags.SOLUTION_VALUE_REQUEST);
                    int requestTimeIndex = inBuf[0], requestDimIndex = inBuf[1];
                    supportPool.execute(() -> {
                        MPI.COMM_WORLD.Isend(new double[]{solution.getValue(requestTimeIndex, requestDimIndex)},
                                0, 1, MPI.DOUBLE, status.source, MpiTags.SOLUTION_VALUE_REPLY);
                    });
                }
            });
            supportThread.setDaemon(true);
            supportThread.start();
            //initial conditions
            parallelSolutionInitialCondition(dimSteps, solution);
        }
        else {
            patternResults = new PatternResult[1];
            patternResults[0] = new PatternResultImpl();
            //solution
            solution = new SolutionArrayImpl(mesh, t0, t1, x0, x1);
            //initial conditions
            sequentialSolutionInitialCondition(dimSteps, solution);
        }

        if (parallelismParam) {
            while (doubleLessOrEquals(t, t1)) {
                if (!(mesh instanceof StaticMesh)) {
                    dimSteps = getDimSteps(t);
                    initMpi(dimSteps.size());
                }

                parallelSolutionSetNextDimRow(timeIndex, t, dimSteps, solution);

                t = mesh.nextTimeStep(t);
                ++timeIndex;
            }

            MPI.COMM_WORLD.Barrier();
            return parallelSolutionCollector.collect(solution, this);
        }
        else {
            while (doubleLessOrEquals(t, t1)) {
                if (!(mesh instanceof StaticMesh)) {
                    dimSteps = getDimSteps(t);
                }

                sequentialSolutionSetNextDimRow(timeIndex, t, dimSteps, solution);

                t = mesh.nextTimeStep(t);
                ++timeIndex;
            }

            return solution;
        }
    }

    protected abstract void sequentialSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps,
                                                            SolutionArrayImpl solution);

    private void sequentialSolutionInitialCondition(List<Double> dimSteps, SolutionArrayImpl solution) {
        double solutionValues[] = new double[dimSteps.size()];

        for (int dimIndex = 0; dimIndex < solutionValues.length; ++dimIndex) {
            solutionValues[dimIndex] = initialDimensionCondition.apply(dimSteps.get(dimIndex));
        }

        solution.setDimRow(solutionValues, 0);
    }

    protected abstract void parallelSolutionSetNextDimRow(int timeIndex, double t, List<Double> dimSteps,
                                                          SolutionArrayImpl solution);

    private void parallelSolutionInitialCondition(List<Double> dimSteps, SolutionArrayImpl solution) {
        double solutionValues[] = new double[counts[rank]];

        try {
            forkJoinPool.submit(() -> IntStream.range(myFromIndex, myToIndex)
                    .parallel()
                    .forEach(dimIndex ->
                            solutionValues[dimIndex - myFromIndex] = initialDimensionCondition.apply(dimSteps.get(dimIndex)))).get();
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }

        solution.setDimRow(solutionValues, 0);
    }

    private boolean doubleLessOrEquals(double a, double b) {
        if (Math.abs(a - b) > 0.0000001) {//a != b
            return a < b;
        }
        else
            return true;
    }

    private List<Double> getDimSteps(double t) {
        double x = x0;
        List<Double> dimSteps = new ArrayList<>();
        while (doubleLessOrEquals(x, x1)) {
            dimSteps.add(x);
            x = mesh.nextDimStep(x, t);
        }
        return dimSteps;
    }

    private void initMpi(int taskSize) {
        int numOfThreads = MPI.COMM_WORLD.Size();
        displacements[0] = 0;
        for (int i = 0; i < numOfThreads; ++i) {
            int iFromIndex = fromIndexOf(i, numOfThreads, taskSize);
            int iToIndex = toIndexOf(i, iFromIndex, numOfThreads, taskSize);
            counts[i] = iToIndex - iFromIndex;
            displacements[i] = iFromIndex;
        }
    }

    private int fromIndexOf(int rank, int numOfThreads, int taskSize) {
        int reminder = taskSize % numOfThreads;
        return rank * (taskSize / numOfThreads) +
                (rank < reminder ? rank : reminder);
    }

    private int toIndexOf(int rank, int fromIndex, int numOfThreads, int taskSize) {
        int reminder = taskSize % numOfThreads;
        return fromIndex + (taskSize / numOfThreads) + (rank < reminder ? 1 : 0);
    }

    protected double[] leftShift(double[] array) {
        double newArray[] = new double[array.length - 1];
        System.arraycopy(array, 1, newArray, 0, newArray.length);
        return newArray;
    }

    protected double[] rightShift(double[] array) {
        double newArray[] = new double[array.length - 1];
        System.arraycopy(array, 0, newArray, 0, newArray.length);
        return newArray;
    }
}
