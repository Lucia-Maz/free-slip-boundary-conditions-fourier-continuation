!======================================================================
! Deformable free surface of order N (Phase 6 of the project that carries
! this patch; verificacion/intent_fase6.txt freezes the conventions).
!
! The exact conditions of a graph surface z = h + eta(x,y,t) are
! transferred to the fixed boundary z = Lz by Taylor series in eta, to
! order N in the amplitude, with the full Navier-Stokes bulk. With
! Tr_M[f] = sum_{m=0}^{M} eta^m/m! d^m f/dz^m (at z = Lz; 0 if M < 0),
! s = grad_h(eta), S = grad v + grad v^T, and p the kinematic pressure
! minus the hydrostatic one at rest:
!
!   kinematic   d(eta)/dt = Tr_{N-1}[w] - s.Tr_{N-2}[u_h]
!   tangential  Tr_{N-1}[S_az] - s_b Tr_{N-2}[S_ab] + s_a Tr_{N-2}[S_zz]
!               - s_a s_b Tr_{N-3}[S_zb] = 0
!   normal      Tr_{N-1}[p] = g eta - gam kappa_N + nu B_N
!     kappa_N = div( s sum_{2j+1<=N} binom(-1/2,j) |s|^{2j} )
!     B_N = sum_{2j<=N-1} (-|s|^2)^j ( Tr_{N-1-2j}[S_zz]
!           - 2 s_a Tr_{N-2-2j}[S_az] + s_a s_b Tr_{N-3-2j}[S_ab] )
!
! These are the Taylor polynomials of the exact conditions of Wang, Tice &
! Kim, ARMA 212 (2014), doi 10.1007/s00205-013-0700-2, sec. 1.1, checked to
! roundoff for N = 1..5 (numerico/fase6/chequeo_expansion.py).
!
! Scheme per Runge-Kutta substage o (fs_general_imposebc): kinematics and
! the explicit part of the normal stress from the previous substage; the
! tangential condition on the new state, its leading term S_az coupled to
! W = w(Lz) exactly (B_k, as in freesurface_imposebc) and its corrections
! as data; the pressure transfer from the pressure being solved. The data
! of the new state are obtained by an affine fixed-point iteration
! accelerated with Anderson mixing.
!
! Normal derivatives of the velocity at the surface are spectral traces
! for every order m. Those of order m >= 2 respond to the Neumann datum
! that is being imposed, as fs_beta(m) per unit datum (measured at setup;
! ~1.40/dz and ~1.15/dz^2 with the default tables), so the tangential
! corrections are coupled to themselves through the pointwise factor
! sigma = beta_2 eta + beta_3 eta^2/2 (N >= 3), and to W through the
! -2 i kx W of the datum. The iteration removes both couplings exactly,
! G = (F - 2 sigma d_x(W_new - W) + sigma Dx)/(1 + sigma), which has the
! same fixed points as F while 1 + sigma > 0. Extrapolating from interior
! nodes instead (no coupling) was tried first and is unstable in time
! at N = 2 for |eta| > ~0.5 dz, and at N = 3 for eta > ~0.4 dz (REGISTRO_fase6.md).
!
! Limit of the order-2 model: where eta < 0 the transferred condition
! d_z u + eta d_zz u = 0 has a spurious boundary-layer solution
! exp((z-Lz)/|eta|), growing as nu/eta^2. On the grid the explicit step
! becomes unstable where 1 + sigma falls to ~0.7 nu dt/dz^2 (eta = -0.667
! dz at nu dt/dz^2 = 0.104; -0.716 dz only as dt -> 0), and the run is
! stopped where 1 + sigma <= 0.05 + nu dt/dz^2. At N = 3 the spurious
! solutions are neutral for uniform eta (weakly growing with slope) and
! 1 + sigma > 0.15; the step is stable except within ~2% of the explicit
! limit (nu dt/dz^2 ~ 0.205, near eta = -1.2 dz), so N = 3 requires
! nu dt/dz^2 <= 0.19.
!
! 2026 L. Mazaira.
!======================================================================

!***********************************************************************
      SUBROUTINE fs_gen_setup(planfc)
!-----------------------------------------------------------------------
!  Surface transforms and trace phases of the general path, and the
!  response fs_beta(m) of the spectral trace d^m u/dz^m at z = Lz to a
!  unit Neumann datum there, measured on a probe field so that it follows
!  the FC-Gram tables in use. Needs only the grids and the FC plan;
!  callable from the tests.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE grid
      USE kes
      USE var

      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: planfc
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: pb
      REAL(KIND=GP) :: xp,ztop
      INTEGER       :: i,j,k,m,p,top

      IF ( ALLOCATED(fs_ex) ) RETURN
      fs_nkx = nx/2
      fs_mx  = 2*nx
      IF ( ny .gt. 1 ) THEN
         fs_my = 2*ny
      ELSE
         fs_my = 1
      ENDIF
      ALLOCATE( fs_ex(fs_mx,fs_nkx), fs_ey(fs_my,ny) )
      DO p = 1,fs_mx
         xp = 2*pi*Lx*real(p-1,kind=GP)/real(fs_mx,kind=GP)
         DO i = 1,fs_nkx
            fs_ex(p,i) = exp(im*kx(i)*xp)
         ENDDO
      ENDDO
      DO p = 1,fs_my
         xp = 2*pi*Ly*real(p-1,kind=GP)/real(fs_my,kind=GP)
         DO j = 1,ny
            fs_ey(p,j) = exp(im*ky(j)*xp)
         ENDDO
      ENDDO

      ! Trace phases at z = Lz, with the 1/nz of the backward FFT. The kz
      ! Nyquist term is the real cosine of the interpolant, whose odd
      ! derivatives vanish at the grid node z = Lz: it is dropped for odd m
      ! (kz = -nz/2 alone would give it a spurious imaginary part).
      top  = nz - Cz
      ztop = z(top)
      ALLOCATE( fs_ph(nz,0:FS_NMAX) )
      DO k = 1,nz
         DO m = 0,FS_NMAX
            fs_ph(k,m) = (im*kz(k))**m*exp(im*kz(k)*ztop)/real(nz,kind=GP)
            IF ( k .eq. nz/2+1 .AND. MOD(m,2) .eq. 1 ) fs_ph(k,m) = 0.0_GP
         ENDDO
      ENDDO

      ! Probe: zero interior, unit Neumann datum at z = Lz
      IF ( .NOT. ALLOCATED(planfc%z%neu) ) THEN
         IF ( myrank .eq. 0 ) THEN
            PRINT*, "[ERROR] fs_gen_setup: the Neumann tables of ", &
                "the z plan are not loaded (load_neumann_tables)."
            FLUSH(6)
         ENDIF
         CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      ENDIF
      ALLOCATE( pb(nz,ny,ista:iend) )
      pb = 0.0_GP
      pb(top,:,:) = 1.0_GP
      CALL neumann_reconstruct(planfc,pb,6,1)
      CALL fftp1d_real_to_complex_z(planfc,pb,MPI_COMM_WORLD)
      DO m = 0,FS_NMAX
         fs_beta(m) = real(SUM(pb(:,1,ista)*fs_ph(:,m)),kind=GP)
      ENDDO
      DEALLOCATE( pb )
      IF ( myrank .eq. 0 ) THEN
         PRINT '(A,4ES12.4)', " [fsorder] trace response to the Neumann datum,"// &
               " beta_m dz^(m-1), m = 0..3:", &
               (fs_beta(m)*(z(2)-z(1))**(m-1), m = 0,FS_NMAX)
      ENDIF

      RETURN
      END SUBROUTINE fs_gen_setup

!***********************************************************************
      SUBROUTINE fs_traces_general(vx,vy,vz,order,tr)
!-----------------------------------------------------------------------
!  Normal derivatives at z = Lz of a velocity field in the (kz,ky,kx)
!  domain, in the (ky,kx) domain with the normalization of (z,ky,kx):
!    tr(:,:,1,m) = d^m u/dz^m, tr(:,:,2,m) = d^m v/dz^m (m = 0..order),
!    tr(:,:,3,m) = d^m w/dz^m.
!  u, v: spectral traces; w: m = 0 by spectral trace, m >= 1 by continuity,
!  d^m w/dz^m = -i (kx d^{m-1}u/dz^{m-1} + ky d^{m-1}v/dz^{m-1}).
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE grid
      USE kes
      USE var

      IMPLICIT NONE

      COMPLEX(KIND=GP),INTENT(IN),DIMENSION(nz,ny,ista:iend) :: vx,vy,vz
      INTEGER, INTENT(IN)  :: order
      COMPLEX(KIND=GP),INTENT(OUT),DIMENSION(ny,ista:iend,3,0:FS_NMAX) :: tr

      COMPLEX(KIND=GP) :: su(0:FS_NMAX),sv(0:FS_NMAX),sw
      INTEGER          :: i,j,k,m

      tr = 0.0_GP
!$omp parallel do if (iend-ista.ge.nth) private (j,k,m,su,sv,sw)
      DO i = ista,iend
         DO j = 1,ny
            su = 0; sv = 0; sw = 0
            DO k = 1,nz
               DO m = 0,order
                  su(m) = su(m) + vx(k,j,i)*fs_ph(k,m)
                  sv(m) = sv(m) + vy(k,j,i)*fs_ph(k,m)
               ENDDO
               sw = sw + vz(k,j,i)*fs_ph(k,0)
            ENDDO
            DO m = 0,order
               tr(j,i,1,m) = su(m)
               tr(j,i,2,m) = sv(m)
            ENDDO
            tr(j,i,3,0) = sw
            DO m = 1,order
               tr(j,i,3,m) = -im*(kx(i)*tr(j,i,1,m-1) + ky(j)*tr(j,i,2,m-1))
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE fs_traces_general

!***********************************************************************
      SUBROUTINE fs_to_phys(nf,Fin,fout)
!-----------------------------------------------------------------------
!  nf surface fields from the (ky,kx) domain (normalization of (z,ky,kx))
!  to the padded physical grid (fs_mx, fs_my), on every rank. Nyquist
!  modes are dropped.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE grid

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nf
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend,nf) :: Fin
      REAL(KIND=GP), INTENT(OUT), DIMENSION(fs_mx,fs_my,nf)    :: fout

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: part,full
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: g
      REAL(KIND=GP) :: fac,ci
      INTEGER :: i,j,l,p,q

      ALLOCATE( part(ny,fs_nkx,nf), full(ny,fs_nkx,nf), g(fs_my,fs_nkx) )
      part = 0.0_GP
      DO l = 1,nf
         DO i = ista,MIN(iend,fs_nkx)
            DO j = 1,ny
               part(j,i,l) = Fin(j,i,l)
            ENDDO
         ENDDO
      ENDDO
      CALL MPI_ALLREDUCE(part,full,ny*fs_nkx*nf,GC_COMPLEX,MPI_SUM, &
                         MPI_COMM_WORLD,ierr)
      IF ( ny .gt. 1 ) full(ny/2+1,:,:) = 0.0_GP
      fac = 1.0_GP/(real(nx,kind=GP)*real(ny,kind=GP))
      DO l = 1,nf
         g = MATMUL(fs_ey,full(:,:,l))
         DO q = 1,fs_my
            DO p = 1,fs_mx
               fout(p,q,l) = real(g(q,1),kind=GP)
            ENDDO
            DO i = 2,fs_nkx
               ci = 2.0_GP
               DO p = 1,fs_mx
                  fout(p,q,l) = fout(p,q,l) + ci*real(fs_ex(p,i)*g(q,i),kind=GP)
               ENDDO
            ENDDO
         ENDDO
         fout(:,:,l) = fout(:,:,l)*fac
      ENDDO
      DEALLOCATE( part,full,g )

      RETURN
      END SUBROUTINE fs_to_phys

!***********************************************************************
      SUBROUTINE fs_to_spec(nf,fin,Fout)
!-----------------------------------------------------------------------
!  Inverse of fs_to_phys: nf padded physical fields to the local (ky,kx)
!  slab, exact projection onto the retained modes (no Nyquist).
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE grid

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nf
      REAL(KIND=GP), INTENT(IN), DIMENSION(fs_mx,fs_my,nf)      :: fin
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend,nf) :: Fout

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: h
      REAL(KIND=GP) :: fac
      INTEGER :: i,j,l,p,q

      ALLOCATE( h(fs_my,fs_nkx) )
      fac = real(nx,kind=GP)*real(ny,kind=GP)/(real(fs_mx,kind=GP)*real(fs_my,kind=GP))
      Fout = 0.0_GP
      DO l = 1,nf
         h = 0.0_GP
         DO i = ista,MIN(iend,fs_nkx)
            DO q = 1,fs_my
               DO p = 1,fs_mx
                  h(q,i) = h(q,i) + fin(p,q,l)*conjg(fs_ex(p,i))
               ENDDO
            ENDDO
         ENDDO
         DO i = ista,MIN(iend,fs_nkx)
            DO j = 1,ny
               IF ( ny .gt. 1 .AND. j .eq. ny/2+1 ) CYCLE
               Fout(j,i,l) = fac*SUM(conjg(fs_ey(:,j))*h(:,i))
            ENDDO
         ENDDO
      ENDDO
      DEALLOCATE( h )

      RETURN
      END SUBROUTINE fs_to_spec

!***********************************************************************
      SUBROUTINE fs_functionals(order,eta,tr,g,gam,nu,K,Dx,Dy,E,Tfull)
!-----------------------------------------------------------------------
!  The order-N surface functionals (see the header of this file), from
!  the elevation eta and the normal derivatives tr of fs_traces_general,
!  all in the (ky,kx) domain:
!    K     kinematic right-hand side
!    Dx,Dy tangential condition minus its leading term S_az(m=0)
!    E     g eta - gam kappa_N + nu B_N (normal stress without pressure)
!    Tfull (optional) the full tangential condition, x and y
!  Products are evaluated on the padded grid and projected exactly.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE grid
      USE kes
      USE var

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: order
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend)             :: eta
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend,3,0:FS_NMAX) :: tr
      REAL(KIND=GP), INTENT(IN) :: g,gam,nu
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend) :: K,Dx,Dy,E
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend,2), OPTIONAL :: Tfull

      ! physical fields: 1 eta, 2 sx, 3 sy, then per m = 0..order a block of
      ! 9: u, v, w, u_x, u_y, v_x, v_y, w_x, w_y
      INTEGER, PARAMETER :: nb = 9
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: Fsp,Gsp
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:)    :: fph,gph
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)      :: et,sx,sy,q2,wk
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:)    :: Sxx,Syy,Sxy,Sxz,Syz,Szz
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:)    :: Uu,Vv,Ww
      REAL(KIND=GP), PARAMETER, DIMENSION(0:3) :: ck = &
          (/ 1.0_GP, -0.5_GP, 0.375_GP, -0.3125_GP /)
      INTEGER :: nf,m,i,j,jj,base,nm

      nm = order
      nf = 3 + nb*(nm+1)
      ALLOCATE( Fsp(ny,ista:iend,nf), fph(fs_mx,fs_my,nf) )
      DO i = ista,iend
         DO j = 1,ny
            Fsp(j,i,1) = eta(j,i)
            Fsp(j,i,2) = im*kx(i)*eta(j,i)
            Fsp(j,i,3) = im*ky(j)*eta(j,i)
            DO m = 0,nm
               base = 3 + nb*m
               Fsp(j,i,base+1) = tr(j,i,1,m)
               Fsp(j,i,base+2) = tr(j,i,2,m)
               Fsp(j,i,base+3) = tr(j,i,3,m)
               Fsp(j,i,base+4) = im*kx(i)*tr(j,i,1,m)
               Fsp(j,i,base+5) = im*ky(j)*tr(j,i,1,m)
               Fsp(j,i,base+6) = im*kx(i)*tr(j,i,2,m)
               Fsp(j,i,base+7) = im*ky(j)*tr(j,i,2,m)
               Fsp(j,i,base+8) = im*kx(i)*tr(j,i,3,m)
               Fsp(j,i,base+9) = im*ky(j)*tr(j,i,3,m)
            ENDDO
         ENDDO
      ENDDO
      CALL fs_to_phys(nf,Fsp,fph)

      ALLOCATE( et(fs_mx,fs_my), sx(fs_mx,fs_my), sy(fs_mx,fs_my) )
      ALLOCATE( q2(fs_mx,fs_my), wk(fs_mx,fs_my) )
      ALLOCATE( Uu(fs_mx,fs_my,0:nm), Vv(fs_mx,fs_my,0:nm), Ww(fs_mx,fs_my,0:nm) )
      ALLOCATE( Sxx(fs_mx,fs_my,0:nm), Syy(fs_mx,fs_my,0:nm), Sxy(fs_mx,fs_my,0:nm) )
      ALLOCATE( Sxz(fs_mx,fs_my,0:nm), Syz(fs_mx,fs_my,0:nm), Szz(fs_mx,fs_my,0:nm) )
      et = fph(:,:,1)
      sx = fph(:,:,2)
      sy = fph(:,:,3)
      q2 = sx**2 + sy**2
      Sxz = 0; Syz = 0; Szz = 0
      DO m = 0,nm
         base = 3 + nb*m
         Uu(:,:,m) = fph(:,:,base+1)
         Vv(:,:,m) = fph(:,:,base+2)
         Ww(:,:,m) = fph(:,:,base+3)
         Sxx(:,:,m) = 2*fph(:,:,base+4)
         Syy(:,:,m) = 2*fph(:,:,base+7)
         Sxy(:,:,m) = fph(:,:,base+5) + fph(:,:,base+6)
      ENDDO
      ! d^m S_az/dz^m = d^{m+1} u_a/dz^{m+1} + d_a d^m w/dz^m ; needs m+1 <= nm
      DO m = 0,nm-1
         base = 3 + nb*m
         Sxz(:,:,m) = Uu(:,:,m+1) + fph(:,:,base+8)
         Syz(:,:,m) = Vv(:,:,m+1) + fph(:,:,base+9)
         Szz(:,:,m) = 2*Ww(:,:,m+1)
      ENDDO

      ALLOCATE( gph(fs_mx,fs_my,6), Gsp(ny,ista:iend,6) )
      ! kinematic
      gph(:,:,1) = trz(Ww,order-1) - sx*trz(Uu,order-2) - sy*trz(Vv,order-2)
      ! tangential, full (x: 2, y: 3), and its leading term (x: 5, y: 6)
      gph(:,:,2) = trz(Sxz,order-1) - (sx*trz(Sxx,order-2) + sy*trz(Sxy,order-2)) &
                + sx*trz(Szz,order-2) - sx*(sx*trz(Sxz,order-3) + sy*trz(Syz,order-3))
      gph(:,:,3) = trz(Syz,order-1) - (sx*trz(Sxy,order-2) + sy*trz(Syy,order-2)) &
                + sy*trz(Szz,order-2) - sy*(sx*trz(Sxz,order-3) + sy*trz(Syz,order-3))
      gph(:,:,5) = Sxz(:,:,0)
      gph(:,:,6) = Syz(:,:,0)
      ! nu B_N
      gph(:,:,4) = 0.0_GP
      jj = 0
      DO WHILE ( 2*jj .le. order-1 )
         wk = trz(Szz,order-1-2*jj) - 2*(sx*trz(Sxz,order-2-2*jj) +  &
              sy*trz(Syz,order-2-2*jj)) + sx*sx*trz(Sxx,order-3-2*jj) + &
              2*sx*sy*trz(Sxy,order-3-2*jj) + sy*sy*trz(Syy,order-3-2*jj)
         gph(:,:,4) = gph(:,:,4) + ((-q2)**jj)*wk
         jj = jj + 1
      ENDDO
      gph(:,:,4) = nu*gph(:,:,4)
      CALL fs_to_spec(6,gph,Gsp)
      K = Gsp(:,:,1)
      Dx = Gsp(:,:,2) - Gsp(:,:,5)
      Dy = Gsp(:,:,3) - Gsp(:,:,6)
      IF ( PRESENT(Tfull) ) THEN
         Tfull(:,:,1) = Gsp(:,:,2)
         Tfull(:,:,2) = Gsp(:,:,3)
      ENDIF

      ! E = g eta - gam kappa_N + nu B_N. kappa_N = lap(eta) for N <= 2 is
      ! taken spectrally; for N >= 3 the flux s sum c_j |s|^2j is projected
      ! and its divergence taken spectrally.
      IF ( order .le. 2 ) THEN
         DO i = ista,iend
            DO j = 1,ny
               E(j,i) = (g + gam*khom(j,i)**2)*eta(j,i) + Gsp(j,i,4)
            ENDDO
         ENDDO
      ELSE
         wk = 0.0_GP
         jj = 0
         DO WHILE ( 2*jj+1 .le. order )
            wk = wk + ck(jj)*q2**jj
            jj = jj + 1
         ENDDO
         gph(:,:,1) = sx*wk
         gph(:,:,2) = sy*wk
         CALL fs_to_spec(2,gph(:,:,1:2),Gsp(:,:,1:2))
         DO i = ista,iend
            DO j = 1,ny
               E(j,i) = g*eta(j,i) - gam*im*(kx(i)*Gsp(j,i,1) + ky(j)*Gsp(j,i,2)) &
                        + Gsp(j,i,4)
            ENDDO
         ENDDO
      ENDIF

      DEALLOCATE( Fsp,fph,gph,Gsp,et,sx,sy,q2,wk,Uu,Vv,Ww,Sxx,Syy,Sxy,Sxz,Syz,Szz )

      CONTAINS

      FUNCTION trz(A,Mt) RESULT(r)
         ! Tr_Mt[A] = sum_{m=0}^{Mt} eta^m/m! A_m (zero if Mt < 0)
         REAL(KIND=GP), INTENT(IN), DIMENSION(fs_mx,fs_my,0:nm) :: A
         INTEGER, INTENT(IN) :: Mt
         REAL(KIND=GP), DIMENSION(fs_mx,fs_my) :: r, pw
         INTEGER :: mm
         r = 0.0_GP
         IF ( Mt .lt. 0 ) RETURN
         pw = 1.0_GP
         DO mm = 0,MIN(Mt,nm)
            r = r + pw*A(:,:,mm)
            pw = pw*et/real(mm+1,kind=GP)
         ENDDO
      END FUNCTION trz

      END SUBROUTINE fs_functionals

!***********************************************************************
      SUBROUTINE fs_ptransfer(order,eta,dder,Pt)
!-----------------------------------------------------------------------
!  Pressure transfer sum_{m=1}^{N-1} eta^m/m! d^m d/dz^m at z = Lz, from
!  the elevation and the normal derivatives dder(:,:,m) of the pressure
!  potential (same units as the potential), in the (ky,kx) domain.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE grid

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: order
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend)           :: eta
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend,0:FS_NMAX) :: dder
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend)          :: Pt

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: Fsp,Gsp
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:)    :: fph,gph
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)      :: pw
      INTEGER :: m

      Pt = 0.0_GP
      IF ( order .lt. 2 ) RETURN
      ALLOCATE( Fsp(ny,ista:iend,order), fph(fs_mx,fs_my,order) )
      ALLOCATE( gph(fs_mx,fs_my,1), Gsp(ny,ista:iend,1), pw(fs_mx,fs_my) )
      Fsp(:,:,1) = eta
      DO m = 1,order-1
         Fsp(:,:,m+1) = dder(:,:,m)
      ENDDO
      CALL fs_to_phys(order,Fsp,fph)
      gph = 0.0_GP
      pw = 1.0_GP
      DO m = 1,order-1
         pw = pw*fph(:,:,1)/real(m,kind=GP)
         gph(:,:,1) = gph(:,:,1) + pw*fph(:,:,m+1)
      ENDDO
      CALL fs_to_spec(1,gph,Gsp)
      Pt = Gsp(:,:,1)
      DEALLOCATE( Fsp,fph,gph,Gsp,pw )

      RETURN
      END SUBROUTINE fs_ptransfer

!***********************************************************************
      SUBROUTINE fs_general_pass(planbc,planfc,vx,vy,vz,pr,dtop0,wst,x,o,ox)
!-----------------------------------------------------------------------
!  One pass of the order-N surface iteration: impose the top datum for
!  the input x = (W, Dx, Dy, Pt), project, and return in ox the outputs
!  (W with the exact B_k update, Dx, Dy from the new state and eta^(o),
!  Pt from the new pressure potential and eta^(o+1)). vx, vy enter in
!  the (z,ky,kx) domain (after the bottom condition) and vz in
!  (kz,ky,kx); they leave in (kz,ky,kx), and pr holds the new potential.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE kes
      USE var

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: vx,vy,vz,pr
      COMPLEX(KIND=GP),INTENT(IN),DIMENSION(planfc%y%n,ista:iend)     :: dtop0,wst
      COMPLEX(KIND=GP),INTENT(IN),DIMENSION(planfc%y%n,ista:iend,4)   :: x
      COMPLEX(KIND=GP),INTENT(OUT),DIMENSION(planfc%y%n,ista:iend,4)  :: ox
      INTEGER, INTENT(IN) :: o

      COMPLEX(KIND=GP), DIMENSION(planfc%y%n,ista:iend) :: dtop,Kd,Ed
      COMPLEX(KIND=GP), DIMENSION(planfc%y%n,ista:iend,0:FS_NMAX) :: dder
      INTEGER :: ind,i,j

      ind = planfc%z%n - planfc%z%c
      DO i = ista,iend
         DO j = 1,planfc%y%n
            vx(ind,j,i) = im*kx(i)*(wst(j,i) - 2*x(j,i,1)) - x(j,i,2)
            vy(ind,j,i) = im*ky(j)*(wst(j,i) - 2*x(j,i,1)) - x(j,i,3)
            dtop(j,i) = dtop0(j,i) - x(j,i,4)
         ENDDO
      ENDDO
      CALL neumann_reconstruct(planfc,vx,6,1)
      CALL neumann_reconstruct(planfc,vy,6,1)
      CALL goto_3d_fourier(planbc,planfc,vx,vy)
      CALL sol_project_fs(vx,vy,vz,pr,dtop,dder)

      ! outputs
      CALL fs_traces_general(vx,vy,vz,fs_order,fs_trn)
      DO i = ista,iend
         DO j = 1,planfc%y%n
            ox(j,i,1) = (fs_trn(j,i,3,0) - fs_bk(j,i)*x(j,i,1))/(1.0_GP - fs_bk(j,i))
         ENDDO
      ENDDO
      CALL fs_functionals(fs_order,fs_eta,fs_trn,fs_grav,fs_tens,fs_nu,Kd, &
                          ox(:,:,2),ox(:,:,3),Ed)
      CALL fs_ptransfer(fs_order,fs_etap,dder,ox(:,:,4))

      RETURN
      END SUBROUTINE fs_general_pass

!***********************************************************************
      SUBROUTINE fs_general_imposebc(planbc,planfc,vx,vy,vz,pr,o,vbot)
!-----------------------------------------------------------------------
!  One Runge-Kutta substage with the order-N deformable surface. See the
!  header of this file. The iteration x -> F(x) over the surface data
!  x = (W, Dx, Dy, Pt) is affine. Its tangential blocks are preconditioned
!  by removing their pointwise self-coupling sigma (header of this file)
!  and the result is accelerated with Anderson mixing (type II, written
!  out in fs_anderson). The blocks are measured against the
!  data they correct: W itself, the full Neumann data of u and v, and the
!  Dirichlet datum of the potential without its mean (the mean is the
!  static g <eta>, trivially satisfied). The iteration stops when input
!  and output agree to fs_tol, or when it has stagnated below 1e2 fs_tol:
!  after at least six passes, the best residual of the last four is not
!  half the best before them. The traces of order m amplify the roundoff
!  of the column as kz^m, which at N = 3 and eta ~ 4 dz leaves a floor of
!  a few 1e-9; the window is long because Anderson is not monotone.
!  The accepted state is the output of the last pass; fs_res reports its
!  residual.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE kes
      USE var
      USE order
      USE grid, ONLY: z

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(IN)  :: planbc
      TYPE(FCPLAN), INTENT(IN)  :: planfc
      COMPLEX(KIND=GP),INTENT(INOUT),DIMENSION(planfc%z%n,planfc%y%n,ista:iend) :: vx,vy,vz,pr
      INTEGER, INTENT(IN)       :: o
      REAL(KIND=GP), INTENT(IN), DIMENSION(2) :: vbot

      COMPLEX(KIND=GP), DIMENSION(planfc%y%n,ista:iend)   :: Kd,Ed,dtop0,wst,wsb,Dx0,Dy0
      COMPLEX(KIND=GP), DIMENSION(planfc%y%n,ista:iend,4) :: x,fx,gx,xn
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: sg
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: ph
      REAL(KIND=GP), DIMENSION(4) :: esc,loc
      REAL(KIND=GP), DIMENSION(fs_maxit) :: rhist
      REAL(KIND=GP) :: fac,res,fac_loc,smin,sthr
      INTEGER       :: i,j,b,it

      IF ( .NOT. fs_ready ) THEN
         IF ( myrank .eq. 0 ) THEN
            PRINT*, "[ERROR] freesurface used before v_fs_init."
            FLUSH(6)
         ENDIF
         CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      ENDIF
      fac = fs_dt/real(o, kind=GP)
      IF ( o .eq. ord ) fs_eta0 = fs_eta

      ! Explicit data from the previous substage (eta^(o+1), traces fs_tr)
      fs_etap = fs_eta
      CALL fs_functionals(fs_order,fs_etap,fs_tr,fs_grav,fs_tens,fs_nu,Kd, &
                          Dx0,Dy0,Ed)
      DO i = ista,iend
         DO j = 1,planfc%y%n
            dtop0(j,i) = fac*Ed(j,i)
            fs_eta(j,i) = fs_eta0(j,i) + fac*Kd(j,i)
         ENDDO
      ENDDO
      ! initial guess of the tangential corrections with eta^(o)
      CALL fs_functionals(fs_order,fs_eta,fs_tr,fs_grav,fs_tens,fs_nu,Kd, &
                          Dx0,Dy0,Ed)

      ! Self-coupling of the tangential corrections through the Neumann
      ! datum, sigma = beta_2 eta + beta_3 eta^2/2, on the padded grid (the
      ! same on every rank)
      ALLOCATE( sg(fs_mx,fs_my), ph(fs_mx,fs_my,6) )
      CALL fs_to_phys(1,fs_eta,ph(:,:,1:1))
      sg = 0.0_GP
      IF ( fs_order .ge. 2 ) sg = fs_beta(2)*ph(:,:,1)
      IF ( fs_order .ge. 3 ) sg = sg + 0.5_GP*fs_beta(3)*ph(:,:,1)**2
      ! guard: at N = 2 the step is unstable where 1 + sigma < ~0.7 nu dt/dz^2
      ! (1D replica, numerico/fase6/estabilidad_cierre.py); margin 0.05
      smin = 1.0_GP + minval(sg)
      sthr = 0.05_GP
      IF ( fs_order .eq. 2 ) sthr = sthr + fs_nu*fs_dt/(z(2)-z(1))**2
      IF ( fs_order .ge. 3 .AND. fs_nu*fs_dt/(z(2)-z(1))**2 .gt. 0.19_GP ) THEN
         IF ( myrank .eq. 0 ) THEN
            PRINT '(A,ES11.3,A)', " [ERROR] freesurface, fsorder = 3 needs"// &
               " nu dt/dz^2 <= 0.19 (got", fs_nu*fs_dt/(z(2)-z(1))**2, "): closer to"// &
               " the explicit limit the spurious modes near eta = -1.2 dz are"// &
               " unstable (fsorder.f90). Reduce dt."
            FLUSH(6)
         ENDIF
         CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      ENDIF
      IF ( smin .le. sthr ) THEN
         IF ( myrank .eq. 0 ) THEN
            PRINT '(A,I2,A,ES11.3,A,ES11.3,A)', &
               " [ERROR] freesurface, fsorder =", fs_order, ": 1 + sigma =", smin, &
               " <= ", sthr, ". The order-2 transfer of the tangential stress"// &
               " is ill-posed where eta < 0 and the step is unstable here"// &
               " (fsorder.f90). Use fsorder = 3, a coarser z grid, or a smaller dt."
            FLUSH(6)
         ENDIF
         CALL MPI_ABORT(MPI_COMM_WORLD,1,ierr)
      ENDIF

      ! Bottom, as in the flat paths; v*_z only through its traces
      CALL fs_trace(planfc,vz,wst)
      IF ( planbc%bczsta .ne. 0 ) CALL fs_trace_bottom(planfc,vz,wsb)
      CALL goto_domain_w_boundaries(planbc,planfc,vx,vy)
      IF ( planbc%bczsta .eq. 0 ) THEN
         CALL noslip_z(planfc,o,vx,vy,vz,pr,vbot,0)
      ELSE
         CALL freeslip_z(planfc,vx,vy,wsb,0)
      ENDIF
      fs_s1 = vx; fs_s2 = vy; fs_s3 = vz

      ! Iteration
      x(:,:,1) = fs_wtop
      x(:,:,2) = Dx0
      x(:,:,3) = Dy0
      x(:,:,4) = fac*fs_ptr_guess
      fs_nhist = 0
      DO it = 1,fs_maxit
         IF ( it .gt. 1 ) THEN
            vx = fs_s1; vy = fs_s2; vz = fs_s3
         ENDIF
         CALL fs_general_pass(planbc,planfc,vx,vy,vz,pr,dtop0,wst,x,o,fx)
         IF ( it .eq. 1 ) THEN
            loc(1) = max(maxval(abs(fx(:,:,1))),maxval(abs(x(:,:,1))))
            loc(2) = 0.0_GP
            DO i = ista,iend
               DO j = 1,planfc%y%n
                  loc(2) = max(loc(2),abs(fx(j,i,2)),abs(fx(j,i,3)),         &
                     abs(im*kx(i)*(wst(j,i)-2*fx(j,i,1))-fx(j,i,2)),        &
                     abs(im*ky(j)*(wst(j,i)-2*fx(j,i,1))-fx(j,i,3)))
               ENDDO
            ENDDO
            loc(3) = loc(2)
            loc(4) = maxval(abs(fx(:,:,4)))
            DO i = ista,iend
               DO j = 1,planfc%y%n
                  IF ( i .eq. 1 .AND. j .eq. 1 ) CYCLE
                  loc(4) = max(loc(4),abs(dtop0(j,i)))
               ENDDO
            ENDDO
            CALL MPI_ALLREDUCE(loc,esc,4,GC_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
            DO b = 1,4
               IF ( esc(b) .le. tiny(1.0_GP) ) esc(b) = 1.0_GP
            ENDDO
         ENDIF
         DO b = 1,4
            loc(b) = maxval(abs(fx(:,:,b)-x(:,:,b)))/esc(b)
         ENDDO
         fac_loc = maxval(loc)
         CALL MPI_ALLREDUCE(fac_loc,res,1,GC_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
         fs_npass = it
         fs_res = res
         rhist(it) = res
         IF ( res .le. fs_tol ) EXIT
         IF ( it .ge. 6 ) THEN
            IF ( res .le. 1.0e2_GP*fs_tol .AND. minval(rhist(it-3:it)) .ge. &
                 0.5_GP*minval(rhist(1:it-4)) ) EXIT
         ENDIF
         ! Preconditioned map on Dx, Dy. The datum of u is i kx (w* - 2W) - Dx
         ! and Dx responds to it as sigma, so with the new W of this pass
         ! (exact, B_k) the fixed point of Dx is
         !   G = (F - 2 sigma d_x(W_new - W) + sigma Dx)/(1 + sigma),
         ! and likewise for Dy with d_y.
         gx = fx
         IF ( fs_order .ge. 2 ) THEN
            DO i = ista,iend
               DO j = 1,planfc%y%n
                  xn(j,i,1) = im*kx(i)*(fx(j,i,1) - x(j,i,1))
                  xn(j,i,4) = im*ky(j)*(fx(j,i,1) - x(j,i,1))
               ENDDO
            ENDDO
            xn(:,:,2:3) = x(:,:,2:3)
            CALL fs_to_phys(4,xn,ph(:,:,1:4))
            CALL fs_to_phys(2,fx(:,:,2:3),ph(:,:,5:6))
            ph(:,:,5) = (ph(:,:,5) - 2*sg*ph(:,:,1) + sg*ph(:,:,2))/(1.0_GP + sg)
            ph(:,:,6) = (ph(:,:,6) - 2*sg*ph(:,:,4) + sg*ph(:,:,3))/(1.0_GP + sg)
            CALL fs_to_spec(2,ph(:,:,5:6),gx(:,:,2:3))
         ENDIF
         CALL fs_anderson(x,gx,esc,xn)
         x = xn
      ENDDO
      DEALLOCATE( sg,ph )
      IF ( it .gt. fs_maxit ) THEN
         fs_nmaxit = fs_nmaxit + 1
         fs_resmaxit = max(fs_resmaxit,fs_res)
      ENDIF

      ! State for the next substage
      fs_wtop = fs_trn(:,:,3,0)
      fs_tr = fs_trn
      fs_ptr_guess = fx(:,:,4)/fac

      RETURN
      END SUBROUTINE fs_general_imposebc

!***********************************************************************
      SUBROUTINE fs_anderson(xa,fxa,esc,xna)
!-----------------------------------------------------------------------
!  One Anderson-mixing update (type II, memory fs_mhist) of the fixed
!  point xa = F(xa) on the surface data, with the blocks scaled by esc:
!  xna = F(xa) - dF gamma, gamma minimizing |r - dR gamma|, r = (F(xa)-xa)/esc.
!  gamma is REAL: the data are Fourier coefficients of real fields and F
!  is real-linear but not complex-linear (products in physical space), so
!  a complex gamma would predict residuals that F does not produce.
!  The history lives in fs_ahx, fs_ahr (module); fs_nhist counts it.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE grid

      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend,4)  :: xa,fxa
      REAL(KIND=GP), INTENT(IN), DIMENSION(4)                  :: esc
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend,4) :: xna

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:,:) :: dR,dF
      COMPLEX(KIND=GP), DIMENSION(ny,ista:iend,4) :: r
      REAL(KIND=GP)    :: Gm(fs_mhist,fs_mhist),bv(fs_mhist),gam(fs_mhist)
      REAL(KIND=GP)    :: loc(fs_mhist*fs_mhist+fs_mhist),glo(fs_mhist*fs_mhist+fs_mhist)
      REAL(KIND=GP) :: tr_
      INTEGER :: b,m,a,c,l,mk

      DO b = 1,4
         r(:,:,b) = (fxa(:,:,b) - xa(:,:,b))/esc(b)
      ENDDO
      ! shift the history: slot 0 is the newest (F, r)
      IF ( fs_nhist .ge. 1 ) THEN
         mk = MIN(fs_nhist,fs_mhist)
         ALLOCATE( dR(ny,ista:iend,4,mk), dF(ny,ista:iend,4,mk) )
         DO m = 1,mk
            IF ( m .eq. 1 ) THEN
               dR(:,:,:,m) = r - fs_ahr(:,:,:,1)
               dF(:,:,:,m) = fxa - fs_ahx(:,:,:,1)
            ELSE
               dR(:,:,:,m) = fs_ahr(:,:,:,m-1) - fs_ahr(:,:,:,m)
               dF(:,:,:,m) = fs_ahx(:,:,:,m-1) - fs_ahx(:,:,:,m)
            ENDIF
         ENDDO
         l = 0
         DO a = 1,mk
            DO c = 1,mk
               l = l + 1
               loc(l) = real(SUM(conjg(dR(:,:,:,a))*dR(:,:,:,c)),kind=GP)
            ENDDO
         ENDDO
         DO a = 1,mk
            l = l + 1
            loc(l) = real(SUM(conjg(dR(:,:,:,a))*r),kind=GP)
         ENDDO
         CALL MPI_ALLREDUCE(loc,glo,l,GC_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
         l = 0
         tr_ = 0.0_GP
         DO a = 1,mk
            DO c = 1,mk
               l = l + 1
               Gm(a,c) = glo(l)
            ENDDO
            tr_ = tr_ + Gm(a,a)
         ENDDO
         DO a = 1,mk
            l = l + 1
            bv(a) = glo(l)
            Gm(a,a) = Gm(a,a) + 1.0e-13_GP*tr_
         ENDDO
         CALL fs_rsolve(mk,Gm(1:mk,1:mk),bv(1:mk),gam(1:mk))
         xna = fxa
         DO m = 1,mk
            xna = xna - gam(m)*dF(:,:,:,m)
         ENDDO
         DEALLOCATE( dR,dF )
      ELSE
         xna = fxa
      ENDIF
      ! push (F, r) into the history
      DO m = fs_mhist,2,-1
         fs_ahx(:,:,:,m) = fs_ahx(:,:,:,m-1)
         fs_ahr(:,:,:,m) = fs_ahr(:,:,:,m-1)
      ENDDO
      fs_ahx(:,:,:,1) = fxa
      fs_ahr(:,:,:,1) = r
      fs_nhist = fs_nhist + 1

      RETURN
      END SUBROUTINE fs_anderson

!***********************************************************************
      SUBROUTINE fs_rsolve(n,A,b,x)
!-----------------------------------------------------------------------
!  Small dense real solve by Gaussian elimination with partial pivoting.
!-----------------------------------------------------------------------
      USE fprecision
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(KIND=GP), INTENT(IN)  :: A(n,n),b(n)
      REAL(KIND=GP), INTENT(OUT) :: x(n)
      REAL(KIND=GP) :: M(n,n),r(n),t,row(n)
      INTEGER :: i,j,p

      M = A
      r = b
      DO j = 1,n
         p = j - 1 + MAXLOC(abs(M(j:n,j)),dim=1)
         IF ( p .ne. j ) THEN
            row = M(j,:); M(j,:) = M(p,:); M(p,:) = row
            t = r(j); r(j) = r(p); r(p) = t
         ENDIF
         IF ( abs(M(j,j)) .le. tiny(1.0_GP) ) THEN
            x = 0.0_GP
            RETURN
         ENDIF
         DO i = j+1,n
            t = M(i,j)/M(j,j)
            M(i,j:n) = M(i,j:n) - t*M(j,j:n)
            r(i) = r(i) - t*r(j)
         ENDDO
      ENDDO
      DO i = n,1,-1
         x(i) = (r(i) - SUM(M(i,i+1:n)*x(i+1:n)))/M(i,i)
      ENDDO

      RETURN
      END SUBROUTINE fs_rsolve

!***********************************************************************
      SUBROUTINE fs_init_pressure(order,eta,tr,pp,ptr)
!-----------------------------------------------------------------------
!  Initial pressure of the order-N surface: p harmonic with dp/dz = 0 at
!  z=0 and p + sum_{m=1}^{N-1} eta^m/m! d^m p/dz^m = E0 at z=Lz, with E0 the
!  normal-stress datum of the initial state, by fixed-point iteration.
!  pp: p in (z,ky,kx); ptr: the converged transfer term (units of p).
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE grid
      USE kes

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: order
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend)             :: eta
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(ny,ista:iend,3,0:FS_NMAX) :: tr
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend)         :: pp
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(ny,ista:iend)            :: ptr

      COMPLEX(KIND=GP), DIMENSION(ny,ista:iend) :: Kd,Dx,Dy,Ed,Pn
      COMPLEX(KIND=GP), DIMENSION(ny,ista:iend,0:FS_NMAX) :: dder
      COMPLEX(KIND=GP), DIMENSION(2,ny,ista:iend) :: bc
      REAL(KIND=GP) :: loc(2),glo(2)
      INTEGER :: i,j,it,top

      top = nz - Cz
      CALL fs_functionals(order,eta,tr,fs_grav,fs_tens,fs_nu,Kd,Dx,Dy,Ed)
      ptr = 0.0_GP
      DO it = 1,60
         DO i = ista,iend
            DO j = 1,ny
               bc(1,j,i) = 0.0_GP
               bc(2,j,i) = Ed(j,i) - ptr(j,i)
            ENDDO
         ENDDO
         CALL laplace_z(bc,pp,fs_work,1,0)
         dder = 0.0_GP
         DO i = ista,iend
            DO j = 1,ny
               dder(j,i,0) = pp(top,j,i)
               dder(j,i,1) = fs_work(top,j,i)
               dder(j,i,2) = khom(j,i)**2*pp(top,j,i)
            ENDDO
         ENDDO
         CALL fs_ptransfer(order,eta,dder,Pn)
         loc(1) = maxval(abs(Pn-ptr))
         loc(2) = maxval(abs(Ed))
         CALL MPI_ALLREDUCE(loc,glo,2,GC_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
         ptr = Pn
         IF ( glo(1) .le. 1.0e-14_GP*max(glo(2),tiny(1.0_GP)) ) EXIT
      ENDDO
      ! final harmonic solve with the converged transfer
      DO i = ista,iend
         DO j = 1,ny
            bc(1,j,i) = 0.0_GP
            bc(2,j,i) = Ed(j,i) - ptr(j,i)
         ENDDO
      ENDDO
      CALL laplace_z(bc,pp,fs_work,1,0)

      RETURN
      END SUBROUTINE fs_init_pressure

!***********************************************************************
      SUBROUTINE fs_order_diagnostic(planfc,a,b,c,t,dt,wallstress)
!-----------------------------------------------------------------------
!  Writes 'freesurface_order_diagnostic.txt' (format frozen in
!  verificacion/intent_fase6.txt):
!    t, eta at k0, 2k0, 3k0 (Re, Im), <eta>, rms(eta),
!    <|T_t|^2>|z=Lz / <|d(v_t)/dz|^2>|z=0 (full order-N tangential
!    condition), passes and final residual of the last substage,
!    max|eta|/dz.
!-----------------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE grid
      USE kes

      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN)  :: planfc
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(IN)    :: dt
      INTEGER, INTENT(IN)          :: t
      DOUBLE PRECISION, INTENT(IN) :: wallstress

      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:,:) :: trd
      COMPLEX(KIND=GP), DIMENSION(ny,ista:iend)   :: Kd,Dx,Dy,Ed
      COMPLEX(KIND=GP), DIMENSION(ny,ista:iend,2) :: Tf
      REAL(KIND=GP), DIMENSION(fs_mx,fs_my,1)     :: ph
      DOUBLE PRECISION :: loc(9),glo(9),fac,wgt,col(13)
      INTEGER :: i,j,n,in,jn,kxn,kyn

      ALLOCATE( trd(ny,ista:iend,3,0:FS_NMAX) )
      CALL fs_traces_general(a,b,c,fs_order,trd)
      CALL fs_functionals(fs_order,fs_eta,trd,fs_grav,fs_tens,fs_nu,Kd,Dx,Dy,Ed,Tf)
      DEALLOCATE( trd )

      loc = 0d0
      DO i = ista,iend
         IF ( i .eq. 1 .OR. i .eq. nx/2+1 ) THEN
            wgt = 1d0
         ELSE
            wgt = 2d0
         ENDIF
         DO j = 1,ny
            loc(8) = loc(8) + wgt*(abs(Tf(j,i,1))**2 + abs(Tf(j,i,2))**2)
            loc(9) = loc(9) + wgt*abs(fs_eta(j,i))**2
         ENDDO
      ENDDO
      IF ( ista .eq. 1 ) loc(7) = real(fs_eta(1,1),kind=GP)
      DO n = 1,3
         kxn = n*fs_kx0
         kyn = n*fs_ky0
         IF ( kxn .gt. nx/2-1 .OR. abs(kyn) .gt. MAX(ny/2-1,0) ) CYCLE
         in = kxn + 1
         IF ( kyn .ge. 0 ) THEN
            jn = kyn + 1
         ELSE
            jn = ny + kyn + 1
         ENDIF
         IF ( in .ge. ista .AND. in .le. iend ) THEN
            loc(2*n-1) = real(fs_eta(jn,in),kind=GP)
            loc(2*n)   = aimag(fs_eta(jn,in))
         ENDIF
      ENDDO
      CALL MPI_REDUCE(loc,glo,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      CALL fs_to_phys(1,fs_eta,ph)
      IF ( myrank .eq. 0 ) THEN
         fac = 1d0/(dble(nx)*dble(ny))
         col(1) = (t-1)*dt
         DO n = 1,3
            IF ( fs_kx0 .eq. 0 .AND. fs_ky0 .eq. 0 ) THEN
               col(2*n) = glo(2*n-1)*fac
               col(2*n+1) = glo(2*n)*fac
            ELSE
               col(2*n) = 2*glo(2*n-1)*fac
               col(2*n+1) = 2*glo(2*n)*fac
            ENDIF
         ENDDO
         col(8) = glo(7)*fac
         col(9) = sqrt(glo(9))*fac
         col(10) = glo(8)*fac**2/max(wallstress,tiny(1d0))
         col(11) = dble(fs_npass)
         col(12) = dble(fs_res)
         col(13) = dble(maxval(abs(ph)))/dble(z(2)-z(1))
         OPEN(1,file='freesurface_order_diagnostic.txt',position='append')
40          FORMAT( 1P 13E23.15 )
            WRITE(1,40) col
         CLOSE(1)
         IF ( fs_nmaxit .gt. 0 ) THEN
            PRINT '(A,I8,A,ES10.2)', " [WARNING] freesurface: ", fs_nmaxit, &
               " substages ended at fsmaxit without converging since the last"// &
               " output; worst final residual", fs_resmaxit
            FLUSH(6)
         ENDIF
      ENDIF
      fs_nmaxit = 0
      fs_resmaxit = 0.0_GP

      RETURN
      END SUBROUTINE fs_order_diagnostic
