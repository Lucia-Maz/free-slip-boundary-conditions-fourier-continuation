! Unit test for the Neumann(z=0) - Dirichlet(z=Lz) branch of laplace_z, the one
! the deformable free surface uses for the pressure.
!
! Solves laplacian(phi) = 0 in a box periodic in x,y with dphi/dz(0) = g0 and
! phi(Lz) = g1, for a handful of (kx,ky) modes, and dumps phi and dphi/dz on the
! physical z grid to 'laplace_neudir.txt'. The check is done outside, by
! verificacion/test_aceptacion_fase5.py, against a closed form in the cosh/sinh
! basis and against the defining properties phi'(0)=g0 and phi(Lz)=g1.
!
! Columns: i, j, khom, z, Re(phi), Im(phi), Re(dphi/dz), Im(dphi/dz),
!          Re(g0), Im(g0), Re(g1), Im(g1)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) THEN
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*- laplace_neudir.f90 *-*-*-*-*-*-*-*-*-*-*-*-*-"
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

C4=0;C5=0;C6=0

! Deterministic complex data, different for every mode and different from the
! ones of laplace_dirneu.f90.
DO i = ista,iend
DO j = 1,ny
   C4(1,j,i) = cmplx(0.7_GP - 0.1_GP*j, 0.3_GP*i, kind=GP)      ! g0 = phi'(0)
   C4(2,j,i) = cmplx(1.2_GP + 0.2_GP*i, -0.4_GP*j, kind=GP)     ! g1 = phi(Lz)
ENDDO
ENDDO

! The (0,0) mode of a real field is real.
IF (ista .eq. 1) THEN
   C4(1,1,1) = cmplx(0.6_GP, 0.0_GP, kind=GP)
   C4(2,1,1) = cmplx(-0.8_GP, 0.0_GP, kind=GP)
ENDIF

CALL laplace_z(C4(1:2,:,:),C5,C6,1,0)

IF (myrank .eq. 0) THEN
   OPEN(1,file='laplace_neudir.txt')
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
   PRINT*, "laplace_neudir: wrote laplace_neudir.txt"
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
