!> @brief Program to analyze a single directed network.
!> @author M. Ángeles Serrano Moral
!>
!> @details This program reads a network from an edge list file, builds an
!>          adjacency list structure, and counts the abundance of various
!>          directed triangle motifs. The final triangle spectrum is written
!>          to an output file.
!
!
! Program to analyze a single directed network.
! Reads an edge list file, builds an adjacency structure, and counts various directed triangle types.
! Writes the triangle spectrum to an output file.

      implicit double precision(x) ! Specifies that variables starting with 'x' are double precision.
      character*80 filename,fout   ! Variables for input and output file names.

      ! Maximum limits for network size and number of edges.
      ! Adjust these parameters if your network exceeds these sizes.
      parameter (NODOSMAX=850000,NEDGESMAX=3600000)
      ! NDEGREEMAXI and NDEGREEMAX2 were removed as they were related to the removed clustering calculation.

      ! Arrays for storing network structure.
      dimension ndegree(1:NODOSMAX)         ! Combined degree (in + out) for each node.
      dimension npunteroini2(1:NODOSMAX)    ! Start index for each node in the combined adjacency list (nvectorlist).
      dimension npunterofin2(1:NODOSMAX)    ! End index (last added element) for each node in the combined list.
      dimension nvectorlist2(1:NEDGESMAX*2)   ! Stores edge direction (+1=out, -1=in, 0=reciprocal). Size needs to accommodate both directions.
      dimension nvectorlist3(1:NEDGESMAX*2)   ! Stores neighbor node ID. Size needs to accommodate both directions.
      dimension internalnet(1:NEDGESMAX,1:2) ! Temporary storage for the input edge list (source, destination).

      ! Initialize degree and pointer arrays to zero.
      data ndegree/NODOSMAX*0/
      data npunteroini2/NODOSMAX*0/
      data npunterofin2/NODOSMAX*0/

      ! Explicitly initialize nvectorlist arrays to zero.
      do i=1,NEDGESMAX*2
        nvectorlist2(i)=0
        nvectorlist3(i)=0
      enddo

      ! Define default input and output file names.
      ! Update filename if your original network file has a different name.
      filename='edge_list.txt'
      fout='edge_list_stat.dat'


      ! Open the input edge list file.
      ! STATUS='UNKNOWN' means create if it doesn't exist, preserve if it does.
      ! Consider using STATUS='OLD' if the file must exist.
      open(1,file=filename,status='unknown')

      NODOS=0   ! Counter for the maximum node ID found (effectively number of nodes if IDs are contiguous starting from 1).
      nlink=0   ! Counter for the number of edges read.

      ! Read edges from the input file.
      do while(.true.)
        ! Read source node (i), destination node (j), and a third value (k, likely weight - unused).
        read(1,*,END=10) i,j,k
        nlink=nlink+1           ! Increment edge count.

        ! Check against maximum edge limit.
        IF (nlink > NEDGESMAX) THEN
            WRITE(*,*) "Error: Exceeded maximum number of edges (NEDGESMAX) during reading."
            STOP
        END IF

        internalnet(nlink,1)=i  ! Store source node ID.
        internalnet(nlink,2)=j  ! Store destination node ID.

        ! Update maximum node ID found.
       if(i.gt.NODOS) NODOS=i
       if(j.gt.NODOS) NODOS=j

      enddo
10    continue ! Label reached upon end-of-file.

      ! Close the input file.
      close(1)

      ! Recalculate degrees after finding the true max node ID.
      ! Initialize degrees to zero before recounting.
      ndegree = 0
      do i=1, nlink
          ! Check against maximum node limit using the final determined NODOS.
          IF (internalnet(i,1) > NODOSMAX .OR. internalnet(i,2) > NODOSMAX) THEN
              WRITE(*,*) "Error: Exceeded maximum number of nodes (NODOSMAX) after reading."
              STOP
          END IF
          ndegree(internalnet(i,1))=ndegree(internalnet(i,1))+1
          ndegree(internalnet(i,2))=ndegree(internalnet(i,2))+1
      enddo


      ! Print basic network statistics to standard output (unit 6).
      write(6,*)'NODOS', NODOS
      write(6,*)'LINKS', nlink

      ! Setup pointers (npunteroini2, npunterofin2) for the combined adjacency list based on degrees.
      ! Pointers cover nodes from 1 to NODOS.
      indexaux=1
       do i=1,NODOS
        if(ndegree(i).gt.0)then ! Only set pointers for nodes with at least one edge.
        npunteroini2(i)=indexaux
        npunterofin2(i)=indexaux-1 ! npunterofin2 points to the last ADDED element (initially before start).
        indexaux=indexaux+ndegree(i)
        ! Check against total degree sum limit (size of nvectorlist*2).
        IF (indexaux-1 > NEDGESMAX*2) THEN
             WRITE(*,*) "Error: Combined degree sum exceeded maximum list capacity (NEDGESMAX*2)."
             STOP
        END IF
        else if (i <= NODOS) then
        ! For isolated nodes within the 1 to NODOS range, set pointers to an empty range.
        npunteroini2(i)=indexaux
        npunterofin2(i)=indexaux-1
        ! write(6,*)'Node with degree 0 ',i ! Diagnostic print for isolated nodes (optional).
        endif
       enddo

       ! Calculate total size used for the combined adjacency list.
       ndegreetotal=indexaux-1
       ! write(6,*)'ndegreetotal',ndegreetotal ! Optional diagnostic print.

       ! Build the combined adjacency list (nvectorlist2, nvectorlist3).
       ! Identifies and marks reciprocal edges (direction 0).
       ! Marks outgoing edges as +1 and incoming edges as -1.
       nreciprocal=0 ! Counter for reciprocal edges (optional, not used in final output).
       do i=1,nlink ! Iterate through each edge from the input file.
          ninside=0
          ! Check if the reverse edge exists in the list for the destination node.
          ! Loop range: from the start pointer to the current end pointer for the destination node.
          ! Ensure the destination node ID is valid.
          IF (internalnet(i,2) > 0 .AND. internalnet(i,2) <= NODOS) THEN
              do j=npunteroini2(internalnet(i,2)), &
     +npunterofin2(internalnet(i,2))
                ! Check index bounds for safety.
                IF (j < 1 .OR. j > NEDGESMAX*2) THEN
                    WRITE(*,*) "Error: Index out of bounds in reciprocal check (j). Index=",j
                    STOP
                END IF
                if(nvectorlist3(j).eq.internalnet(i,1))then
                  ninside=1        ! Reciprocal edge found.
                  nvectorlist2(j)=0 ! Mark the edge in destination's list as reciprocal.
                  nreciprocal=nreciprocal+1 ! Increment reciprocal counter.
                  ! Find the original edge in the source's list and mark it as reciprocal.
                  ! Loop range: from the start pointer to the current end pointer for the source node.
                  ! Ensure the source node ID is valid.
                  IF (internalnet(i,1) > 0 .AND. internalnet(i,1) <= NODOS) THEN
                      do l=npunteroini2(internalnet(i,1)), &
           +npunterofin2(internalnet(i,1))
                         ! Check index bounds for safety.
                         IF (l < 1 .OR. l > NEDGESMAX*2) THEN
                             WRITE(*,*) "Error: Index out of bounds in reciprocal check (l). Index=",l
                             STOP
                         END IF
                        if(nvectorlist3(l).eq.internalnet(i,2))then
                          nvectorlist2(l)=0 ! Mark the edge in source's list as reciprocal.
                          goto 30           ! Branch after finding/marking reciprocal.
                        endif
                      enddo
                  END IF
                endif
              enddo
          END IF
30        continue ! Label for GOTO.

          ! If reciprocal edge was not found, add both directed edges to the lists.
          if(ninside.eq.0)then
            ! Add outgoing edge from source.
            ! Ensure source node ID is valid.
            IF (internalnet(i,1) > 0 .AND. internalnet(i,1) <= NODOS) THEN
                npunterofin2(internalnet(i,1))= &
     +               npunterofin2(internalnet(i,1))+1
                ! Check index bounds before writing.
                IF (npunterofin2(internalnet(i,1)) < 1 .OR. npunterofin2(internalnet(i,1)) > NEDGESMAX*2) THEN
                     WRITE(*,*) "Error: Index out of bounds adding outgoing edge. Index=",npunterofin2(internalnet(i,1))
                     STOP
                END IF
                nvectorlist3(npunterofin2(internalnet(i,1)))= &
     +               internalnet(i,2)
                nvectorlist2(npunterofin2(internalnet(i,1)))=1 ! Direction +1 (outgoing).
            END IF

             ! Add incoming edge to destination (represented as outgoing from dest to source).
             ! Ensure destination node ID is valid.
             IF (internalnet(i,2) > 0 .AND. internalnet(i,2) <= NODOS) THEN
                 npunterofin2(internalnet(i,2))= &
      +               npunterofin2(internalnet(i,2))+1
                 ! Check index bounds before writing.
                 IF (npunterofin2(internalnet(i,2)) < 1 .OR. npunterofin2(internalnet(i,2)) > NEDGESMAX*2) THEN
                     WRITE(*,*) "Error: Index out of bounds adding incoming edge. Index=",npunterofin2(internalnet(i,2))
                     STOP
                 END IF
                 nvectorlist3(npunterofin2(internalnet(i,2)))= &
      +               internalnet(i,1)
                 nvectorlist2(npunterofin2(internalnet(i,2)))=-1 ! Direction -1 (incoming).
             END IF
           endif
       enddo


!%%%%%%%%%%%%%%%%%%%%%%%%%Triangle Counting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Open the output file for writing triangle counts.
      open(2,file=fout,status='unknown')

      ! Initialize triangle type counters (accumulate triplets*3 as in original).
      nc1=0 ! 3-cycle
      nc2=0 ! 3-nocycle
      nc3=0 ! 4-1biflow
      nc4=0 ! 4-1biout (relative to the bidirectional edge in patterns)
      nc5=0 ! 4-1biin (relative to the bidirectional edge in patterns)
      nc6=0 ! 5-2bi (Triangle with two reciprocal edges, one directed)
      nc7=0 ! 6-3bi (Triangle with three reciprocal edges)
      nctot=0 ! Total triplets found (each triangle counted 3 times).

      ! Main loop to iterate through nodes and count triangles (triplets).
      do i=1,NODOS
      ! Node 'i' can be the center of a triplet only if it has degree > 1.
      if(ndegree(i).gt.1)then
        nclust=0 ! Reset triplet counter for node 'i'.
        ! Iterate over distinct pairs of neighbors j and k of node i using list indices.
        do j_idx=npunteroini2(i),npunterofin2(i)-1
           ! Check index bounds for safety.
           IF (j_idx < 1 .OR. j_idx > NEDGESMAX*2) THEN
              WRITE(*,*) "Error: Index out of bounds in triplet search (j_idx). Index=",j_idx
              STOP
           END IF
           do k_idx=j_idx+1,npunterofin2(i)
              ! Check index bounds for safety.
              IF (k_idx < 1 .OR. k_idx > NEDGESMAX*2) THEN
                 WRITE(*,*) "Error: Index out of bounds in triplet search (k_idx). Index=",k_idx
                 STOP
              END IF

               ! Get the IDs of the neighbor nodes of 'i'.
               neighbor_j = nvectorlist3(j_idx)
               neighbor_k = nvectorlist3(k_idx)

               ! Search for the third edge between neighbor j and neighbor k.
               ! Iterate over the neighbors of neighbor j.
               ! Check if neighbor node ID is valid before accessing its pointers.
               IF (neighbor_j > 0 .AND. neighbor_j <= NODOS) THEN
                   do l_idx=npunteroini2(neighbor_j), &
        +              npunterofin2(neighbor_j)
                       ! Check index bounds for safety.
                       IF (l_idx < 1 .OR. l_idx > NEDGESMAX*2) THEN
                           WRITE(*,*) "Error: Index out of bounds in triplet search (l_idx). Index=",l_idx
                           STOP
                       END IF

                       ! If neighbor k is found in neighbor j's list, a triplet (i, neighbor_j, neighbor_k) is found.
                       if(nvectorlist3(l_idx).eq.neighbor_k)then
                         nclust=nclust+1 ! Increment triplet count for node 'i'.
                         ! Get directionality flags for the three edges in the triplet.
                         ! nt1: direction of (i, neighbor_j) from i's list.
                         ! nt2: direction of (i, neighbor_k) from i's list.
                         ! nt3: direction of (neighbor_j, neighbor_k) from neighbor_j's list.
                         nt1=nvectorlist2(j_idx)
                         nt2=nvectorlist2(k_idx)
                         nt3=nvectorlist2(l_idx)
                         ! Classify the triangle based on edge directions (ported logic).
                         if(abs(nt1)+abs(nt2).eq.0)then ! (i,j) and (i,k) are reciprocal.
                           if(nt3.eq.0)then ! (j,k) is reciprocal.
                             nc7=nc7+1 ! 6-3bi (All three edges reciprocal).
                           else ! (j,k) is directed.
                             nc6=nc6+1 ! 5-2bi (Two reciprocal, one directed).
                           endif

                         elseif(abs(nt1)+abs(nt2).eq.1)then ! One reciprocal, one directed from i.
                           if(nt3.eq.0)then ! (j,k) is reciprocal.
                             nc6=nc6+1 ! 5-2bi.
                           else ! (j,k) is directed.
                             if(nt2.eq.0)then ! (i,k) reciprocal, (i,j) directed.
                               if(nt1*nt3.eq.1)then
                                 nc3=nc3+1 ! 4-1biflow (like i->j->k<->i).
                               else
                                 if(nt1.eq.1)then
                                   nc4=nc4+1 ! 4-1biout.
                                 else
                                   nc5=nc5+1 ! 4-1biin.
                                 endif
                               endif
                             else ! nt1 .eq. 0: (i,j) reciprocal, (i,k) directed.
                               if(nt2*nt3.eq.-1)then
                                 nc3=nc3+1 ! 4-1biflow (like i<->j->k<-i).
                               else
                                 if(nt2.eq.1)then
                                   nc4=nc4+1 ! 4-1biout.
                                 else
                                   nc5=nc5+1 ! 4-1biin.
                                 endif
                               endif
                             endif
                           endif

                         elseif(abs(nt1)+abs(nt2).eq.2)then ! Both (i,j) and (i,k) are directed from i.
                           if(nt3.eq.0)then ! (j,k) is reciprocal.
                             if(nt1+nt2.eq.2)then ! i->j, i->k.
                               nc5=nc5+1 ! 4-1biin (Two outgoing from i, reciprocal base).
                             elseif(nt1+nt2.eq.-2)then ! i<-j, i<-k.
                               nc4=nc4+1 ! 4-1biout (Two incoming to i, reciprocal base).
                             else ! i->j, i<-k OR i<-j, i->k.
                               nc3=nc3+1 ! 4-1biflow.
                             endif
                           else ! (j,k) is directed.
                             if((nt1+nt2.eq.2).or.(nt1+nt2.eq.-2))then ! i->j, i->k OR i<-j, i<-k.
                               nc2=nc2+1 ! 3-nocycle.
                             else ! i->j, i<-k OR i<-j, i->k.
                               if(nt1*nt3.eq.1)then
                                 nc1=nc1+1 ! 3-cycle.
                               else
                                 nc2=nc2+1 ! 3-nocycle.
                               endif
                             endif
                           endif

                         endif ! End of classification IFs.
                       endif ! End IF (third link found).
                   enddo ! End DO l_idx (search in neighbor_j's list).
               END IF ! End IF (neighbor node ID is valid).
           enddo ! End DO k_idx (second neighbor of i).
        enddo ! End DO j_idx (first neighbor of i).


        ! Accumulate total number of triplets found (each triangle counted 3 times).
        nctot=nctot+nclust

      endif ! End IF (node has degree > 1).
      enddo ! End DO i (loop over all nodes).

!%%%%%%%%%%%%%%%%%%%%%End Triangle Counting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! Write final triangle counts (triplets / 3) to the output file (unit 2).
       write(2,*)'clustering spectrum'
       write(2,*)'3-cycle ', INT(dble(nc1)/3.0)
       write(2,*)'3-nocycle ', INT(dble(nc2)/3.0)
       write(2,*)'4-1biflow ', INT(dble(nc3)/3.0)
       write(2,*)'4-1biout ', INT(dble(nc4)/3.0)
       write(2,*)'4-1biin ', INT(dble(nc5)/3.0)
       write(2,*)'5-2bi ', INT(dble(nc6)/3.0)
       write(2,*)'6-3bi ', INT(dble(nc7)/3.0)
       write(2,*)'total number of triangles ', INT(dble(nctot)/3.0)
      close(2) ! Close the output file.


      stop ! Program termination.

      end ! End of the main program.