! Unit test for the Dirichlet(z=0) - Neumann(z=Lz) branch of laplace_z.
!
! Solves laplacian(phi) = 0 in a box periodic in x,y with phi(0) = g0 and
! dphi/dz(Lz) = g1, for a handful of (kx,ky) modes, and dumps phi and dphi/dz on
! the physical z grid to 'laplace_dirneu.txt'. The check itself is done outside,
! by verificacion/test_aceptacion_fase2.py, against a closed form written in the
! cosh/sinh basis -- a different algebraic route from the exp(-2 k Lz) one the code
! uses -- and against the defining properties phi(0)=g0 and phi'(Lz)=g1.
!
! Columns: i, j, khom, z, Re(phi), Im(phi), Re(dphi/dz), Im(dphi/dz),
!          Re(g0), Im(g0), Re(g1), Im(g1)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) THEN
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*- laplace_dirneu.f90 *-*-*-*-*-*-*-*-*-*-*-*-*-"
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

C4=0;C5=0;C6=0

! Boundary data. Deterministic, complex, and different for every mode so that a
! branch that silently ignored one of the two data would be caught.
DO i = ista,iend
DO j = 1,ny
   C4(1,j,i) = cmplx(1.0_GP + 0.1_GP*i, 0.2_GP*j, kind=GP)     ! g0 = phi(0)
   C4(2,j,i) = cmplx(-0.3_GP*j, 0.5_GP + 0.1_GP*i, kind=GP)    ! g1 = phi'(Lz)
ENDDO
ENDDO

! The (kx,ky) = (0,0) mode of a real field is real, and laplace_z takes the real
! part there on purpose (every branch does). Feeding it complex data would be
! unphysical, and the branch would legitimately drop the imaginary part.
IF (ista .eq. 1) THEN
   C4(1,1,1) = cmplx(1.1_GP, 0.0_GP, kind=GP)
   C4(2,1,1) = cmplx(-0.3_GP, 0.0_GP, kind=GP)
ENDIF

CALL laplace_z(C4(1:2,:,:),C5,C6,0,1)

IF (myrank .eq. 0) THEN
   OPEN(1,file='laplace_dirneu.txt')
   DO i = ista,MIN(iend,3)
   DO j = 1,MIN(ny,3)
   DO k = 1,nz-Cz
      WRITE(1,*) i, j, khom(j,i), z(k), &
                 real(C5(k,j,i), kind=GP), aimag(C5(k,j,i)), &
                 real(C6(k,j,i), kind=GP), aimag(C6(k,j,i)), &
                 real(C4(1,j,i), kind=GP), aimag(C4(1,j,i)), &
                 real(C4(2,j,i), kind=GP), aimag(C4(2,j,i))
   ENDDO
   ENDDO
   ENDDO
   CLOSE(1)
   PRINT*, "laplace_dirneu: wrote laplace_dirneu.txt"
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
