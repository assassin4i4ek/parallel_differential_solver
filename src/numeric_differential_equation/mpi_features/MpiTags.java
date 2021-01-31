package numeric_differential_equation.mpi_features;

public class MpiTags {
    public static final int RIGHT_COLUMN_FIRST_ELEMENT = 0;
    public static final int RIGHT_COLUMN_LAST_ELEMENT = 1;
    public static final int MAIN_DIAGONAL_FIRST_ELEMENT = 2;
    public static final int MAIN_DIAGONAL_LAST_ELEMENT = 3;
    public static final int LEFT_COLUMN_FIRST_ELEMENT = 4;
    public static final int LEFT_COLUMN_LAST_ELEMENT = 5;
    public static final int FIRST_FREE_TERM = 6;
    public static final int LAST_FREE_TERM = 7;
    public static final int UPPER_BOUND_RESULT = 8;
    public static final int LOWER_BOUND_RESULT = 9;
    public static final int SOLUTION_VALUE_REQUEST = 10;
    public static final int SOLUTION_VALUE_REPLY = 11;
    public static final int ARRAY_APPROXIMATELY_ZERO = 12;
    private MpiTags() {}
}
