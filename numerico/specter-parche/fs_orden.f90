! Unit tests of the order-N deformable surface (Phase 6): the normal
! derivatives at z = Lz and the order-N surface functionals, for analytic
! solenoidal fields frozen in verificacion/intent_fase6.txt. The checks are
! done outside, by verificacion/test_aceptacion_fase6.py (H2, H3), against
! the exact geometry evaluated at the displaced surface.
!
!   fs_trazas.txt : m, i, j, Re, Im of d^m u/dz^m at z = Lz (m = 0..3)
!   fs_jet.txt    : N, eps, id, i, j, Re, Im, id = 1 K, 2 T_x, 3 T_y,
!                   4 Phi = Tr_{N-1}[p] - g eta + gam kappa_N - nu B_N,
!                   5 kappa_N

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF (myrank .eq. 0) THEN
PRINT*, "-*-*-*-*-*-*-*-*-*-*-*-*-*- fs_orden.f90 *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-"
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

BLOCK
   COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:,:) :: trj,trz0
   COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:)   :: pdj,Tfj
   COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:)     :: etaj,Kj,Dxj,Dyj,Ej,Ptj,Fj
   REAL(KIND=GP) :: epsj(4) = (/ 0.2_GP, 0.1_GP, 0.05_GP, 0.025_GP /)
   REAL(KIND=GP) :: ampj
   INTEGER :: nj,ie,idf,ii,jj,kk,mm

   ! the tests skip setup_bc, which is where the main binary loads the
   ! Neumann tables of the z plan; fs_gen_setup needs them for its probe
   IF ( .NOT. ALLOCATED(planfc%z%neu) ) CALL load_neumann_tables(planfc%z,planfc%tdir,1)
   CALL fs_gen_setup(planfc)
   ALLOCATE( trj(ny,ista:iend,3,0:FS_NMAX), trz0(ny,ista:iend,3,0:FS_NMAX) )
   ALLOCATE( pdj(ny,ista:iend,0:FS_NMAX), Tfj(ny,ista:iend,2) )
   ALLOCATE( etaj(ny,ista:iend), Kj(ny,ista:iend), Dxj(ny,ista:iend) )
   ALLOCATE( Dyj(ny,ista:iend), Ej(ny,ista:iend), Ptj(ny,ista:iend), Fj(ny,ista:iend) )

   ! Fields in the physical domain
   DO kk = ksta,pkend
   DO jj = 1,ny
   DO ii = 1,nx
      R1(ii,jj,kk) = -0.3_GP*sin(2*x(ii)+y(jj))*(1+0.3_GP*z(kk)**2) &
                     + 0.48_GP*sin(x(ii)+y(jj))*sin(1.2_GP*z(kk))
      R2(ii,jj,kk) = 0.3_GP*cos(x(ii)-y(jj))*exp(0.6_GP*z(kk)) &
                     + 0.6_GP*sin(2*x(ii)+y(jj))*(1+0.3_GP*z(kk)**2)
      R3(ii,jj,kk) = 0.4_GP*cos(x(ii)+y(jj))*cos(1.2_GP*z(kk)) &
                     - 0.5_GP*sin(x(ii)-y(jj))*exp(0.6_GP*z(kk))
   ENDDO
   ENDDO
   ENDDO
   CALL fftp3d_real_to_complex(planfc,R1,C1,MPI_COMM_WORLD)
   CALL fftp3d_real_to_complex(planfc,R2,C2,MPI_COMM_WORLD)
   CALL fftp3d_real_to_complex(planfc,R3,C3,MPI_COMM_WORLD)
   DO kk = ksta,pkend
   DO jj = 1,ny
   DO ii = 1,nx
      R1(ii,jj,kk) = 0.6_GP*cos(x(ii))*cosh(1.1_GP*z(kk)) &
                     + 0.3_GP*sin(y(jj)-x(ii))*exp(0.4_GP*z(kk))
   ENDDO
   ENDDO
   ENDDO
   CALL fftp3d_real_to_complex(planfc,R1,C4,MPI_COMM_WORLD)

   ! Normal derivatives
   CALL fs_traces_general(C1,C2,C3,FS_NMAX,trj)
   pdj = 0.0_GP
   DO ii = ista,iend
   DO jj = 1,ny
   DO mm = 0,2
      pdj(jj,ii,mm) = SUM(C4(:,jj,ii)*fs_ph(:,mm))
   ENDDO
   ENDDO
   ENDDO
   IF (myrank .eq. 0) THEN
      OPEN(1,file='fs_trazas.txt')
      DO mm = 0,3
      DO ii = ista,MIN(iend,3)
      DO jj = 1,MIN(ny,3)
         WRITE(1,'(3I5,2ES26.16)') mm, ii, jj, real(trj(jj,ii,1,mm),kind=GP), &
                                   aimag(trj(jj,ii,1,mm))
      ENDDO
      ENDDO
      ENDDO
      CLOSE(1)
   ENDIF

   ! Order-N functionals for eta = eps [cos(x) + 0.5 sin(x+2y)]
   trz0 = 0.0_GP
   IF (myrank .eq. 0) OPEN(2,file='fs_jet.txt')
   DO nj = 1,3
   DO ie = 1,4
      etaj = 0.0_GP
      ampj = epsj(ie)*real(nx,kind=GP)*real(ny,kind=GP)
      IF ( 2 .ge. ista .AND. 2 .le. iend ) THEN
         etaj(1,2) = etaj(1,2) + ampj/2
         etaj(3,2) = etaj(3,2) - im*0.25_GP*ampj
      ENDIF
      CALL fs_functionals(nj,etaj,trj,1.3_GP,0.7_GP,0.11_GP,Kj,Dxj,Dyj,Ej,Tfj)
      CALL fs_ptransfer(nj,etaj,pdj,Ptj)
      DO idf = 1,5
         SELECT CASE (idf)
         CASE (1)
            Fj = Kj
         CASE (2)
            Fj = Tfj(:,:,1)
         CASE (3)
            Fj = Tfj(:,:,2)
         CASE (4)
            Fj = pdj(:,:,0) + Ptj - Ej
         CASE (5)
            CALL fs_functionals(nj,etaj,trz0,0.0_GP,1.0_GP,0.0_GP,Kj,Dxj,Dyj,Fj)
            Fj = -Fj
         END SELECT
         IF (myrank .eq. 0) THEN
            DO ii = ista,MIN(iend,nx/2)
            DO jj = 1,ny
               IF ( ny .gt. 1 .AND. jj .eq. ny/2+1 ) CYCLE
               WRITE(2,'(I3,ES14.6,3I5,2ES26.16)') nj, epsj(ie), idf, ii, jj, &
                     real(Fj(jj,ii),kind=GP), aimag(Fj(jj,ii))
            ENDDO
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   ENDDO
   IF (myrank .eq. 0) THEN
      CLOSE(2)
      PRINT*, "fs_orden: wrote fs_trazas.txt and fs_jet.txt"
   ENDIF
   DEALLOCATE( trj,trz0,pdj,Tfj,etaj,Kj,Dxj,Dyj,Ej,Ptj,Fj )
END BLOCK
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
