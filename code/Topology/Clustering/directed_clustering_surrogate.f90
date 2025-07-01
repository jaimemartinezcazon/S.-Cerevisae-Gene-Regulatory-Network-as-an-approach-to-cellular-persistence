!> @brief Program to analyze an ensemble of surrogate directed networks.
!> @author Jaime Martínez Cazón
!>
!> @details This program iterates through a collection of edge list files, each
!>          representing a null model network. For each network, it counts the
!>          abundance of directed triangle motifs. Finally, it calculates the
!>          mean and standard deviation for each motif type across the entire
!>          ensemble and writes the statistics to an output file.
!
!
! Analyzes edge list files generated for null models.
! Calculates means and standard deviations of triangle counts.
! Incorporates triangle counting logic ported from original user Fortran code.

PROGRAM analyze_surrogates

  IMPLICIT NONE ! Explicitly declare all variables

  ! --- Configuration Parameters ---
  ! Directory containing the .txt edge list files for null models.
  ! Update this path according to your file location.
  CHARACTER(LEN=200), PARAMETER :: EDGE_LIST_DIR = "Null_Models_EdgeLists"
  ! Total number of null models to process.
  ! Ensure this matches the number of files in the specified directory.
  INTEGER, PARAMETER :: NUM_NULL_MODELS = 1000
  ! Name of the output file for statistics.
  ! Change if a different output file name is desired.
  CHARACTER(LEN=200), PARAMETER :: OUTPUT_FILE = "surrogate_list_stat.dat"

  ! Maximum estimated size for temporary edge arrays when reading a single .txt file.
  ! Increase this value if your null networks have significantly more edges.
  ! An "out of bounds" error during initial reading indicates this value is too small.
  ! Dynamic arrays are used for the network structure within the counting subroutine.
  INTEGER, PARAMETER :: MAX_EDGES_READ = 3700000

  ! Arrays to store triangle counts for each null model.
  ! Dimensioned by the total number of models. Element 'i' corresponds to model 'i'.
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_3cycle
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_3nocycle
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_4_1biflow
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_4_1biout
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_4_1biin
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_5_2bi
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_6_3bi ! Corresponds to nc7 in original logic
  INTEGER, DIMENSION(NUM_NULL_MODELS) :: count_total_triangles ! Sum of the above for each model (divided by 3)

  ! Variables to store calculated means and standard deviations.
  REAL :: mean_3cycle, stddev_3cycle
  REAL :: mean_3nocycle, stddev_3nocycle
  REAL :: mean_4_1biflow, stddev_4_1biflow
  REAL :: mean_4_1biout, stddev_4_1biout
  REAL :: mean_4_1biin, stddev_4_1biin
  REAL :: mean_5_2bi, stddev_5_2bi
  REAL :: mean_6_3bi, stddev_6_3bi
  REAL :: mean_total_triangles, stddev_total_triangles

  ! Loop index, file unit numbers, and I/O status variables.
  INTEGER :: i                ! Loop counter for null models.
  INTEGER :: file_unit        ! File unit number for reading .txt files.
  INTEGER :: ioerr            ! Variable to capture I/O errors (e.g., file open).
  INTEGER :: read_stat        ! Variable to capture read statement status.
  CHARACTER(LEN=250) :: filename ! String to construct the full file path.

  ! Temporary arrays to hold edge data read from one file.
  ! These are re-populated for each null model.
  INTEGER, DIMENSION(MAX_EDGES_READ) :: src_nodes    ! Source node IDs (1-based).
  INTEGER, DIMENSION(MAX_EDGES_READ) :: dest_nodes   ! Destination node IDs (1-based).
  REAL, DIMENSION(MAX_EDGES_READ) :: weights        ! Edge weights (read but not used in counting).
  INTEGER :: num_edges_this_graph             ! Actual number of edges read from the current file.

  ! Open the output file for writing statistics.
  ! STATUS='REPLACE' ensures the file is created or overwritten.
  OPEN(UNIT=10, FILE=TRIM(OUTPUT_FILE), STATUS='REPLACE', ACTION='WRITE', IOSTAT=ioerr)
  IF (ioerr /= 0) THEN
    WRITE(*,*) "Error opening output file:", TRIM(OUTPUT_FILE)
    STOP ! Terminate if output file cannot be opened.
  END IF

  WRITE(*,*) "Starting analysis of", NUM_NULL_MODELS, "null models..."
  WRITE(*,*) "Reading files from:", TRIM(EDGE_LIST_DIR)

  ! --- Process Each Null Model ---
  file_unit = 20 ! Assign a file unit for reading input files.

  DO i = 1, NUM_NULL_MODELS
    ! Construct the full file path for the current model.
    ! Uses I4.4 format for zero-padding model numbers (e.g., 0001, 0002).
    WRITE(filename, '(A,A,I4.4,A)') TRIM(EDGE_LIST_DIR), '/null_model_', i, '.txt'

    ! Open the current model's edge list file.
    OPEN(UNIT=file_unit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', &
         FORM='FORMATTED', IOSTAT=ioerr)

    IF (ioerr /= 0) THEN
      ! Report error and skip if the file cannot be opened.
      WRITE(*,*) "Warning: Could not open file", TRIM(filename)
      WRITE(*,*) "Counters for this model will be set to 0."
      ! Initialize counters for the skipped model to zero.
      count_3cycle(i) = 0
      count_3nocycle(i) = 0
      count_4_1biflow(i) = 0
      count_4_1biout(i) = 0
      count_4_1biin(i) = 0
      count_5_2bi(i) = 0
      count_6_3bi(i) = 0
      count_total_triangles(i) = 0
      ! No need to close a file if OPEN failed with STATUS='OLD'.
      CYCLE ! Proceed to the next model.
    END IF

    ! --- Read Edges from File ---
    num_edges_this_graph = 0

    DO
        ! Read one edge (source, target, weight).
        ! Branch to label 200 upon reaching end-of-file.
        READ(file_unit, *, IOSTAT=read_stat, END=200) src_nodes(num_edges_this_graph+1), &
                                                      dest_nodes(num_edges_this_graph+1), &
                                                      weights(num_edges_this_graph+1)

        ! Check for read errors other than end-of-file.
        IF (read_stat < 0) THEN
            WRITE(*,*) "Unexpected read error in file", TRIM(filename)
            WRITE(*,*) "read_stat =", read_stat, "at approximate edge", num_edges_this_graph + 1
            CLOSE(file_unit)
            CLOSE(10)
            STOP
        END IF

        ! Increment edge counter upon successful read.
        num_edges_this_graph = num_edges_this_graph + 1

        ! Check for array overflow during reading.
        IF (num_edges_this_graph > MAX_EDGES_READ) THEN
           WRITE(*,*) "Error: Exceeded MAX_EDGES_READ limit."
           WRITE(*,*) "num_edges_this_graph =", num_edges_this_graph
           WRITE(*,*) "MAX_EDGES_READ =", MAX_EDGES_READ
           WRITE(*,*) "Increase MAX_EDGES_READ parameter and recompile."
           CLOSE(file_unit)
           CLOSE(10)
           STOP
        END IF

    END DO ! End of read loop

200 CONTINUE ! Label for end-of-file branch.

    CLOSE(file_unit) ! Close the current input file.

    ! num_edges_this_graph now holds the total number of edges read.

    ! --- Count Triangles for the Current Network ---
    ! Call the subroutine to perform triangle counting based on the edge list.
    ! Results for this model are stored in the i-th element of the count arrays.
    CALL count_triangles_sub(src_nodes, dest_nodes, num_edges_this_graph, &
                             count_3cycle(i), count_3nocycle(i), &
                             count_4_1biflow(i), count_4_1biout(i), count_4_1biin(i), &
                             count_5_2bi(i), count_6_3bi(i), count_total_triangles(i))


    ! Print progress periodically.
    IF (MOD(i, 100) == 0) THEN
      WRITE(*,*) "Processed model", i, "of", NUM_NULL_MODELS
    END IF

  END DO ! End of loop over all null models

  WRITE(*,*) "Analysis of", NUM_NULL_MODELS, "models completed."
  WRITE(*,*) "Calculating statistics (Mean and Standard Deviation)..."

  ! --- Calculate Statistics (Mean and Standard Deviation) ---
  ! Calculate the mean for each triangle type count.
  mean_3cycle = calculate_mean(count_3cycle, NUM_NULL_MODELS)
  mean_3nocycle = calculate_mean(count_3nocycle, NUM_NULL_MODELS)
  mean_4_1biflow = calculate_mean(count_4_1biflow, NUM_NULL_MODELS)
  mean_4_1biout = calculate_mean(count_4_1biout, NUM_NULL_MODELS)
  mean_4_1biin = calculate_mean(count_4_1biin, NUM_NULL_MODELS)
  mean_5_2bi = calculate_mean(count_5_2bi, NUM_NULL_MODELS)
  mean_6_3bi = calculate_mean(count_6_3bi, NUM_NULL_MODELS)
  mean_total_triangles = calculate_mean(count_total_triangles, NUM_NULL_MODELS)

  ! Calculate the sample standard deviation for each triangle type count.
  ! Requires at least two data points (NUM_NULL_MODELS > 1).
  IF (NUM_NULL_MODELS > 1) THEN
    stddev_3cycle = calculate_stddev(count_3cycle, mean_3cycle, NUM_NULL_MODELS)
    stddev_3nocycle = calculate_stddev(count_3nocycle, mean_3nocycle, NUM_NULL_MODELS)
    stddev_4_1biflow = calculate_stddev(count_4_1biflow, mean_4_1biflow, NUM_NULL_MODELS)
    stddev_4_1biout = calculate_stddev(count_4_1biout, mean_4_1biout, NUM_NULL_MODELS)
    stddev_4_1biin = calculate_stddev(count_4_1biin, mean_4_1biin, NUM_NULL_MODELS)
    stddev_5_2bi = calculate_stddev(count_5_2bi, mean_5_2bi, NUM_NULL_MODELS)
    stddev_6_3bi = calculate_stddev(count_6_3bi, mean_6_3bi, NUM_NULL_MODELS)
    stddev_total_triangles = calculate_stddev(count_total_triangles, mean_total_triangles, NUM_NULL_MODELS)
  ELSE
    ! Standard deviation is not applicable for N=1.
    stddev_3cycle = 0.0
    stddev_3nocycle = 0.0
    stddev_4_1biflow = 0.0
    stddev_4_1biout = 0.0
    stddev_4_1biin = 0.0
    stddev_5_2bi = 0.0
    stddev_6_3bi = 0.0
    stddev_total_triangles = 0.0
    WRITE(*,*) "Only 1 null model processed. Standard deviation is 0."
  END IF

  ! --- Write Results to Output File ---
  WRITE(10,*) "clustering spectrum"
  ! Write formatted output matching the example layout.
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "3-cycle           ", INT(mean_3cycle), "+/-", stddev_3cycle
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "3-nocycle         ", INT(mean_3nocycle), "+/-", stddev_3nocycle
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "4-1biflow         ", INT(mean_4_1biflow), "+/-", stddev_4_1biflow
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "4-1biout          ", INT(mean_4_1biout), "+/-", stddev_4_1biout
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "4-1biin           ", INT(mean_4_1biin), "+/-", stddev_4_1biin
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "5-2bi             ", INT(mean_5_2bi), "+/-", stddev_5_2bi
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "6-3bi             ", INT(mean_6_3bi), "+/-", stddev_6_3bi
  WRITE(10,*) "" ! Blank line for separation.
  WRITE(10, '(A, 1X, I10, 1X, A, 1X, F8.3)') "total number of triangles", INT(mean_total_triangles), "+/-", stddev_total_triangles

  CLOSE(10) ! Close the output file.

  WRITE(*,*) "Results saved to", TRIM(OUTPUT_FILE)
  WRITE(*,*) "Program finished."

CONTAINS

  ! --- Function to Calculate Mean of an Integer Array ---
  FUNCTION calculate_mean(data, n) RESULT(mean_val)
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: data ! Input array of integers.
    INTEGER, INTENT(IN) :: n                   ! Number of elements in the array.
    REAL :: mean_val                           ! Result: the mean value (as a real).
    INTEGER :: i
    REAL :: sum_val

    sum_val = 0.0
    ! Sum all elements, converting to real.
    DO i = 1, n
      sum_val = sum_val + REAL(data(i))
    END DO

    ! Calculate the mean. Handle division by zero if n is 0.
    IF (n > 0) THEN
      mean_val = sum_val / REAL(n)
    ELSE
      mean_val = 0.0
    END IF
  END FUNCTION calculate_mean

  ! --- Function to Calculate Sample Standard Deviation ---
  FUNCTION calculate_stddev(data, mean_val, n) RESULT(stddev_val)
    IMPLICIT NONE
    INTEGER, DIMENSION(:), INTENT(IN) :: data     ! Input array of integers.
    REAL, INTENT(IN) :: mean_val                  ! The pre-calculated mean of the data.
    INTEGER, INTENT(IN) :: n                      ! Number of elements.
    REAL :: stddev_val                            ! Result: the standard deviation.
    INTEGER :: i
    REAL :: sum_sq_diff                           ! Sum of squared differences from the mean.

    ! Sample standard deviation requires at least 2 data points (N > 1).
    IF (n <= 1) THEN
      stddev_val = 0.0 ! Standard deviation is undefined or 0.
      RETURN           ! Exit function.
    END IF

    sum_sq_diff = 0.0
    ! Sum the squared differences between each data point and the mean.
    DO i = 1, n
      sum_sq_diff = sum_sq_diff + (REAL(data(i)) - mean_val)**2
    END DO

    ! Calculate standard deviation: square root of the variance (sum of squared diffs / (N-1)).
    stddev_val = SQRT(sum_sq_diff / REAL(n - 1))

  END FUNCTION calculate_stddev


  ! --- Subroutine for Triangle Counting for a Single Network ---
  ! Ports the triangle counting logic from the original code.
  ! Takes the edge list for the current network and returns the 8 triangle counts.

  SUBROUTINE count_triangles_sub(src_in, dest_in, num_edges_in, &
                                 count_3cycle_out, count_3nocycle_out, &
                                 count_4_1biflow_out, count_4_1biout_out, count_4_1biin_out, &
                                 count_5_2bi_out, count_6_3bi_out, count_total_out)

    IMPLICIT NONE ! Explicitly declare all local variables

    ! Input Arguments: Edge data for THIS network.
    INTEGER, DIMENSION(:), INTENT(IN) :: src_in    ! Array of source node IDs (1-based).
    INTEGER, DIMENSION(:), INTENT(IN) :: dest_in   ! Array of destination node IDs (1-based).
    ! REAL, DIMENSION(:), INTENT(IN) :: wt_in     ! Array of weights (not used in counting logic).
    INTEGER, INTENT(IN) :: num_edges_in          ! Number of edges in the input arrays.

    ! Output Arguments: Calculated triangle counts for THIS network.
    INTEGER, INTENT(OUT) :: count_3cycle_out      ! Count for 3-cycle type.
    INTEGER, INTENT(OUT) :: count_3nocycle_out    ! Count for 3-nocycle type.
    INTEGER, INTENT(OUT) :: count_4_1biflow_out   ! Count for 4-1biflow type.
    INTEGER, INTENT(OUT) :: count_4_1biout_out    ! Count for 4-1biout type.
    INTEGER, INTENT(OUT) :: count_4_1biin_out     ! Count for 4-1biin type.
    INTEGER, INTENT(OUT) :: count_5_2bi_out       ! Count for 5-2bi type.
    INTEGER, INTENT(OUT) :: count_6_3bi_out       ! Count for 6-3bi type (Based on original nc7).
    INTEGER, INTENT(OUT) :: count_total_out       ! Sum of the above counts (divided by 3).

    ! --- Local Variables (Ported from Original Code) ---
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ndegree_local ! Combined degree (in+out) for each node.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: npunteroini2_local ! Start index in nvectorlist for each node.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: npunterofin2_local ! End index (last added element) in nvectorlist.
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nvectorlist2_local ! Direction (+1=out, -1=in, 0=reciprocal).
    INTEGER, ALLOCATABLE, DIMENSION(:) :: nvectorlist3_local ! Neighbor node ID.

    INTEGER :: current_NODOS      ! Number of nodes in this network.
    INTEGER :: local_nlink        ! Number of edges (same as num_edges_in).
    INTEGER :: max_node_id        ! Highest node ID found (used for array sizing).
    INTEGER :: i, j, k, l         ! Loop indices.
    INTEGER :: j_idx, k_idx, l_idx ! Indices within nvectorlist.
    INTEGER :: neighbor_j, neighbor_k ! IDs of neighbor nodes.

    INTEGER :: indexaux           ! Auxiliary index for building pointers.
    INTEGER :: ndegreetotal_local ! Total size required for nvectorlist (sum of all degrees).
    INTEGER :: ninside            ! Flag to check for reciprocity during list construction.
    INTEGER :: nreciprocal_local  ! Counter for reciprocal edges (not used in final counts, but ported).
    INTEGER :: nclust             ! Counter for triplets centered at node 'i'.
    INTEGER :: nt1, nt2, nt3      ! Directionality flags for edges in a triplet.

    ! Local counters for triangle types (accumulate triplets * 3, like original nc1-nc7).
    INTEGER :: local_nc1 ! 3-cycle
    INTEGER :: local_nc2 ! 3-nocycle
    INTEGER :: local_nc3 ! 4-1biflow
    INTEGER :: local_nc4 ! 4-1biout
    INTEGER :: local_nc5 ! 4-1biin
    INTEGER :: local_nc6 ! 5-2bi
    INTEGER :: local_nc7 ! 6-3bi (Based on original nc7)
    INTEGER :: local_nctot ! Total found triplets * 3.

    INTEGER :: ioerr_alloc ! Status for allocation errors.

    ! --- Step 1: Determine Number of Nodes and Calculate Degrees ---
    local_nlink = num_edges_in
    max_node_id = 0

    ! First pass to find the maximum node ID.
    DO i = 1, local_nlink
      max_node_id = MAX(max_node_id, src_in(i), dest_in(i))
    END DO
    current_NODOS = max_node_id

    ! Handle case of empty or single-node network (no triangles).
    IF (current_NODOS <= 1) THEN
        count_3cycle_out = 0
        count_3nocycle_out = 0
        count_4_1biflow_out = 0
        count_4_1biout_out = 0
        count_4_1biin_out = 0
        count_5_2bi_out = 0
        count_6_3bi_out = 0
        count_total_out = 0
        RETURN ! Exit subroutine.
    END IF

    ! Allocate memory for degree array.
    ALLOCATE(ndegree_local(1:current_NODOS), STAT=ioerr_alloc)
    IF (ioerr_alloc /= 0) THEN
        WRITE(*,*) "Error allocating ndegree_local, stat=", ioerr_alloc
        STOP ! Terminate on allocation failure.
    END IF
    ndegree_local = 0 ! Initialize degrees to zero.

    ! Calculate combined degrees (in+out).
    DO i = 1, local_nlink
      ! Check node IDs are within the allocated range.
      IF (src_in(i) < 1 .OR. src_in(i) > current_NODOS .OR. &
          dest_in(i) < 1 .OR. dest_in(i) > current_NODOS) THEN
          WRITE(*,*) "Error: Node IDs out of range in input data."
          WRITE(*,*) "Edge:", src_in(i), dest_in(i), " Max ID:", current_NODOS
          IF (ALLOCATED(ndegree_local)) DEALLOCATE(ndegree_local)
          STOP
      END IF
      ndegree_local(src_in(i)) = ndegree_local(src_in(i)) + 1
      ndegree_local(dest_in(i)) = ndegree_local(dest_in(i)) + 1
    END DO

    ! Calculate total degree sum, needed for nvectorlist size.
    ndegreetotal_local = 0
    DO i = 1, current_NODOS
        ndegreetotal_local = ndegreetotal_local + ndegree_local(i)
    END DO

    ! Handle case of no edges after finding max node ID.
    IF (ndegreetotal_local <= 0) THEN
       count_3cycle_out = 0
       count_3nocycle_out = 0
       count_4_1biflow_out = 0
       count_4_1biout_out = 0
       count_4_1biin_out = 0
       count_5_2bi_out = 0
       count_6_3bi_out = 0
       count_total_out = 0
       IF (ALLOCATED(ndegree_local)) DEALLOCATE(ndegree_local)
       RETURN
    END IF

    ! Allocate memory for pointer arrays and combined lists.
    ALLOCATE(npunteroini2_local(1:current_NODOS), STAT=ioerr_alloc)
    IF (ioerr_alloc /= 0) THEN
        WRITE(*,*) "Error allocating npunteroini2_local, stat=", ioerr_alloc
        STOP
    END IF
    ALLOCATE(npunterofin2_local(1:current_NODOS), STAT=ioerr_alloc)
    IF (ioerr_alloc /= 0) THEN
        WRITE(*,*) "Error allocating npunterofin2_local, stat=", ioerr_alloc
        STOP
    END IF
    ALLOCATE(nvectorlist2_local(1:ndegreetotal_local), STAT=ioerr_alloc)
    IF (ioerr_alloc /= 0) THEN
        WRITE(*,*) "Error allocating nvectorlist2_local, stat=", ioerr_alloc
        STOP
    END IF
    ALLOCATE(nvectorlist3_local(1:ndegreetotal_local), STAT=ioerr_alloc)
    IF (ioerr_alloc /= 0) THEN
        WRITE(*,*) "Error allocating nvectorlist3_local, stat=", ioerr_alloc
        STOP
    END IF


    ! Initialize pointer and list arrays.
    npunteroini2_local = 0
    npunterofin2_local = 0
    nvectorlist2_local = 0
    nvectorlist3_local = 0

    ! --- Step 2: Build Combined Adjacency Structure with Directionality ---
    ! Populate pointer arrays.
    indexaux = 1
    DO i = 1, current_NODOS
        npunteroini2_local(i) = indexaux
        npunterofin2_local(i) = indexaux - 1 ! npunterofin2 points to the last ADDED element
        indexaux = indexaux + ndegree_local(i)
    END DO

    nreciprocal_local = 0
    ! Iterate over each original edge (source, dest).
    DO i = 1, local_nlink
        ! Check if the reverse edge (dest_in(i), src_in(i)) already exists in the combined list
        ! being built for dest_in(i).
        ninside = 0
        DO j = npunteroini2_local(dest_in(i)), npunterofin2_local(dest_in(i))
            IF (j < 1 .OR. j > ndegreetotal_local) THEN
              WRITE(*,*) "Index error j in reciprocal search:", j, " total=", ndegreetotal_local
              GOTO 999 ! Clean up and stop
            END IF
            IF (nvectorlist3_local(j) .eq. src_in(i)) THEN
                ninside = 1 ! Reciprocal edge found.
                ! Mark both directions as reciprocal (0 in nvectorlist2).
                nvectorlist2_local(j) = 0 ! Edge (dest_in(i), src_in(i)) in dest_in(i)'s list.
                nreciprocal_local = nreciprocal_local + 1

                ! Find the original edge (src_in(i), dest_in(i)) in src_in(i)'s list
                ! and mark it as reciprocal.
                DO l = npunteroini2_local(src_in(i)), npunterofin2_local(src_in(i))
                   IF (l < 1 .OR. l > ndegreetotal_local) THEN
                      WRITE(*,*) "Index error l in reciprocal search:", l, " total=", ndegreetotal_local
                      GOTO 999 ! Clean up and stop
                   END IF
                    IF (nvectorlist3_local(l) .eq. dest_in(i)) THEN
                        nvectorlist2_local(l) = 0 ! Edge (src_in(i), dest_in(i)) in src_in(i)'s list.
                        GOTO 30 ! Branch after finding/marking reciprocal edge.
                    END IF
                END DO ! End DO l
                GOTO 30 ! Branch if reverse found but forward not (shouldn't happen with correct logic).
            END IF
        END DO ! End DO j (search in dest_in(i)'s list)
30        CONTINUE ! Label for GOTO.


        ! If reciprocal edge was NOT found, add both directions to lists.
        IF (ninside .eq. 0) THEN
            ! Add (src_in(i), dest_in(i)) to src_in(i)'s list.
            npunterofin2_local(src_in(i)) = npunterofin2_local(src_in(i)) + 1
            IF (npunterofin2_local(src_in(i)) > ndegreetotal_local) THEN
               WRITE(*,*) "Index error adding outgoing edge (src, dest). Position exceeded:", &
                          npunterofin2_local(src_in(i)), " Total allocated:", ndegreetotal_local
               GOTO 999 ! Clean up and stop
            END IF
            nvectorlist3_local(npunterofin2_local(src_in(i))) = dest_in(i)
            nvectorlist2_local(npunterofin2_local(src_in(i))) = 1 ! Outgoing from source.

            ! Add (dest_in(i), src_in(i)) to dest_in(i)'s list (as an "inverted" edge).
            npunterofin2_local(dest_in(i)) = npunterofin2_local(dest_in(i)) + 1
             IF (npunterofin2_local(dest_in(i)) > ndegreetotal_local) THEN
               WRITE(*,*) "Index error adding incoming edge (dest, src). Position exceeded:", &
                          npunterofin2_local(dest_in(i)), " Total allocated:", ndegreetotal_local
               GOTO 999 ! Clean up and stop
             END IF
            nvectorlist3_local(npunterofin2_local(dest_in(i))) = src_in(i)
            nvectorlist2_local(npunterofin2_local(dest_in(i))) = -1 ! Incoming to destination.
        END IF
    END DO ! End DO i (loop over original edges)


    ! --- Step 3: Count Triangles Using the Constructed Structure ---
    ! Ported directly from the original counting logic.
    local_nc1 = 0 ! 3-cycle
    local_nc2 = 0 ! 3-nocycle
    local_nc3 = 0 ! 4-1biflow
    local_nc4 = 0 ! 4-1biout
    local_nc5 = 0 ! 4-1biin
    local_nc6 = 0 ! 5-2bi
    local_nc7 = 0 ! 6-3bi (Based on original nc7)
    local_nctot = 0 ! Total found triplets (each triangle counted 3 times)

    DO i = 1, current_NODOS
        ! A node 'i' can be the center of a triplet only if it has at least two entries
        ! in its combined adjacency list.
        IF ((npunterofin2_local(i) - npunteroini2_local(i) + 1) .ge. 2) THEN
            nclust = 0 ! Reset triplet counter for node 'i'.

            ! Iterate over distinct pairs of neighbors of 'i' (neighbor_j, neighbor_k)
            ! using indices j_idx and k_idx in nvectorlist.
            DO j_idx = npunteroini2_local(i), npunterofin2_local(i) - 1
               IF (j_idx < 1 .OR. j_idx > ndegreetotal_local) THEN
                   WRITE(*,*) "Index error j_idx loop, j_idx=",j_idx
                   GOTO 999 ! Clean up and stop
               END IF

               DO k_idx = j_idx + 1, npunterofin2_local(i)
                 IF (k_idx < 1 .OR. k_idx > ndegreetotal_local) THEN
                   WRITE(*,*) "Index error k_idx loop, k_idx=",k_idx
                   GOTO 999 ! Clean up and stop
                 END IF

                   ! Get the IDs of the neighbor nodes of 'i'.
                   neighbor_j = nvectorlist3_local(j_idx)
                   neighbor_k = nvectorlist3_local(k_idx)

                   ! Check if node neighbor_k is a neighbor of node neighbor_j.
                   ! Search for neighbor_k in the combined list of neighbor_j.
                   ! Ensure neighbor_j is a valid node ID and has entries in its list.
                   IF (neighbor_j > 0 .AND. neighbor_j <= current_NODOS .AND. &
                       (npunterofin2_local(neighbor_j) - npunteroini2_local(neighbor_j) + 1) .ge. 1) THEN

                       ! Iterate over the neighbors of neighbor_j.
                       DO l_idx = npunteroini2_local(neighbor_j), npunterofin2_local(neighbor_j)
                           IF (l_idx < 1 .OR. l_idx > ndegreetotal_local) THEN
                               WRITE(*,*) "Index error l_idx loop, l_idx=",l_idx
                               GOTO 999 ! Clean up and stop
                           END IF

                           ! If neighbor_k is found in neighbor_j's list, a triplet (i, neighbor_j, neighbor_k) is found.
                           IF (nvectorlist3_local(l_idx) .eq. neighbor_k) THEN
                               nclust = nclust + 1 ! Increment triplet count centered at i.

                               ! Get directionality flags for the three edges in the triplet.
                               ! nt1: directionality of (i, neighbor_j) from i's list.
                               ! nt2: directionality of (i, neighbor_k) from i's list.
                               ! nt3: directionality of (neighbor_j, neighbor_k) from neighbor_j's list.
                               nt1 = nvectorlist2_local(j_idx)
                               nt2 = nvectorlist2_local(k_idx)
                               nt3 = nvectorlist2_local(l_idx)

                               ! Classify the triangle based on directionalities.
                               ! Ported directly from the original code's IF/ELSEIF structure.
                               IF(abs(nt1)+abs(nt2) .eq. 0)THEN ! (i,j) and (i,k) are reciprocal.
                                 IF(nt3 .eq. 0)THEN ! (j,k) is reciprocal.
                                   local_nc7=local_nc7+1 ! 6-3bi (All three edges reciprocal).
                                 ELSE ! (j,k) is directed.
                                   local_nc6=local_nc6+1 ! 5-2bi (Two reciprocal, one directed).
                                 ENDIF

                               ELSEIF(abs(nt1)+abs(nt2) .eq. 1)THEN ! One reciprocal, one directed from i.
                                 IF(nt3 .eq. 0)THEN ! (j,k) is reciprocal.
                                   local_nc6=local_nc6+1 ! 5-2bi.
                                 ELSE ! (j,k) is directed.
                                   IF(nt2 .eq. 0)THEN ! (i,k) reciprocal, (i,j) directed.
                                     IF(nt1*nt3 .eq. 1)THEN
                                       local_nc3=local_nc3+1 ! 4-1biflow (like i->j->k<->i).
                                     ELSE
                                       IF(nt1 .eq. 1)THEN
                                         local_nc4=local_nc4+1 ! 4-1biout.
                                       ELSE
                                         local_nc5=local_nc5+1 ! 4-1biin.
                                       ENDIF
                                     ENDIF
                                   ELSE ! nt1 .eq. 0: (i,j) reciprocal, (i,k) directed.
                                     IF(nt2*nt3 .eq. -1)THEN
                                       local_nc3=local_nc3+1 ! 4-1biflow (like i<->j->k<-i).
                                     ELSE
                                       IF(nt2 .eq. 1)THEN
                                         local_nc4=local_nc4+1 ! 4-1biout.
                                       ELSE
                                         local_nc5=local_nc5+1 ! 4-1biin.
                                       ENDIF
                                     ENDIF
                                   ENDIF
                                 ENDIF

                               ELSEIF(abs(nt1)+abs(nt2) .eq. 2)THEN ! Both (i,j) and (i,k) are directed from i.
                                 IF(nt3 .eq. 0)THEN ! (j,k) is reciprocal.
                                   IF(nt1+nt2 .eq. 2)THEN ! i->j, i->k.
                                     local_nc5=local_nc5+1 ! 4-1biin (Two outgoing, one reciprocal base).
                                   ELSEIF(nt1+nt2 .eq. -2)THEN ! i<-j, i<-k.
                                     local_nc4=local_nc4+1 ! 4-1biout (Two incoming, one reciprocal base).
                                   ELSE ! i->j, i<-k OR i<-j, i->k.
                                     local_nc3=local_nc3+1 ! 4-1biflow.
                                   ENDIF
                                 ELSE ! (j,k) is directed.
                                   IF((nt1+nt2 .eq. 2) .OR. (nt1+nt2 .eq. -2))THEN ! i->j, i->k OR i<-j, i<-k.
                                     local_nc2=local_nc2+1 ! 3-nocycle.
                                   ELSE ! i->j, i<-k OR i<-j, i->k.
                                     IF(nt1*nt3 .eq. 1)THEN
                                       local_nc1=local_nc1+1 ! 3-cycle.
                                     ELSE
                                       local_nc2=local_nc2+1 ! 3-nocycle.
                                     ENDIF
                                   ENDIF
                                 ENDIF

                               ENDIF ! End of classification IFs.

                           ENDIF ! End IF (third link found).
                       ENDDO ! End DO l_idx (search in neighbor_j's list).
                   ENDIF ! End IF (neighbor_j is valid and has neighbors).
               ENDDO ! End DO k_idx (second neighbor of i).
            ENDDO ! End DO j_idx (first neighbor of i).

            ! Accumulate triplets centered at i to total triplets.
            local_nctot = local_nctot + nclust

        ENDIF ! End IF (node i has enough neighbors).
    ENDDO ! End DO i (loop over all nodes).

    ! The counters local_nc1...local_nc7 and local_nctot represent triplets * 3.
    ! Divide by 3.0 to get the actual triangle counts, as in the original output format.

    ! --- Step 4: Assign Results to Output Variables ---
    count_3cycle_out = INT(dble(local_nc1) / 3.0)
    count_3nocycle_out = INT(dble(local_nc2) / 3.0)
    count_4_1biflow_out = INT(dble(local_nc3) / 3.0)
    count_4_1biout_out = INT(dble(local_nc4) / 3.0)
    count_4_1biin_out = INT(dble(local_nc5) / 3.0)
    count_5_2bi_out = INT(dble(local_nc6) / 3.0)
    count_6_3bi_out = INT(dble(local_nc7) / 3.0)
    count_total_out = INT(dble(local_nctot) / 3.0)


    ! --- Deallocate Dynamically Allocated Memory ---
999 CONTINUE ! Label for error cleanup branch.
    ! Deallocate only if the arrays were allocated.
    IF (ALLOCATED(ndegree_local)) THEN
        DEALLOCATE(ndegree_local, STAT=ioerr_alloc)
        IF (ioerr_alloc /= 0) THEN
            WRITE(*,*) "Error deallocating ndegree_local, stat=", ioerr_alloc
        END IF
    END IF

    IF (ALLOCATED(npunteroini2_local)) THEN
        DEALLOCATE(npunteroini2_local, STAT=ioerr_alloc)
        IF (ioerr_alloc /= 0) THEN
            WRITE(*,*) "Error deallocating npunteroini2_local, stat=", ioerr_alloc
        END IF
    END IF

    IF (ALLOCATED(npunterofin2_local)) THEN
        DEALLOCATE(npunterofin2_local, STAT=ioerr_alloc)
        IF (ioerr_alloc /= 0) THEN
            WRITE(*,*) "Error deallocating npunterofin2_local, stat=", ioerr_alloc
        END IF
    END IF

    IF (ALLOCATED(nvectorlist2_local)) THEN
        DEALLOCATE(nvectorlist2_local, STAT=ioerr_alloc)
        IF (ioerr_alloc /= 0) THEN
            WRITE(*,*) "Error deallocating nvectorlist2_local, stat=", ioerr_alloc
        END IF
    END IF

    IF (ALLOCATED(nvectorlist3_local)) THEN
        DEALLOCATE(nvectorlist3_local, STAT=ioerr_alloc)
        IF (ioerr_alloc /= 0) THEN
            WRITE(*,*) "Error deallocating nvectorlist3_local, stat=", ioerr_alloc
        END IF
    END IF

  END SUBROUTINE count_triangles_sub

END PROGRAM analyze_surrogates