MODULE FORCING_FIELDS_TYPE_MOD

   USE PARKIND_WAVE, ONLY : JWRB, JWIM, JWRO
   USE FIELD_MODULE, ONLY : FIELD_3RB, FIELD_2IM, FIELD_2RB, FIELD_3IM, FIELD_4RB, FIELD_4IM, FIELD_2RD, &
   &                        FIELD_2RM, FIELD_1IM, FIELD_1RB
   USE FIELD_FACTORY_MODULE, ONLY : FIELD_NEW, FIELD_DELETE
   IMPLICIT NONE

   PRIVATE

   INTERFACE ALLOC
      MODULE PROCEDURE FORCING_FIELDS_ALLOC
   END INTERFACE ALLOC
   INTERFACE DEALLOC
      MODULE PROCEDURE FORCING_FIELDS_DEALLOC
   END INTERFACE DEALLOC
   INTERFACE SYNC_DEVICE_RDWR
      MODULE PROCEDURE FORCING_FIELDS_SYNC_DEVICE_RDWR
   END INTERFACE SYNC_DEVICE_RDWR
   INTERFACE SYNC_DEVICE_RDONLY
      MODULE PROCEDURE FORCING_FIELDS_SYNC_DEVICE_RDONLY
   END INTERFACE SYNC_DEVICE_RDONLY
   INTERFACE SYNC_HOST_RDWR
      MODULE PROCEDURE FORCING_FIELDS_SYNC_HOST_RDWR
   END INTERFACE SYNC_HOST_RDWR
   INTERFACE SYNC_HOST_RDONLY
      MODULE PROCEDURE FORCING_FIELDS_SYNC_HOST_RDONLY
   END INTERFACE SYNC_HOST_RDONLY
   INTERFACE GET_DEVICE_DATA_RDWR
      MODULE PROCEDURE FORCING_FIELDS_GET_DEVICE_DATA_RDWR
   END INTERFACE GET_DEVICE_DATA_RDWR
   INTERFACE GET_DEVICE_DATA_RDONLY
      MODULE PROCEDURE FORCING_FIELDS_GET_DEVICE_DATA_RDONLY
   END INTERFACE GET_DEVICE_DATA_RDONLY
   INTERFACE GET_DEVICE_DATA_WRONLY
      MODULE PROCEDURE FORCING_FIELDS_GET_DEVICE_DATA_WRONLY
   END INTERFACE GET_DEVICE_DATA_WRONLY
   INTERFACE GET_HOST_DATA_RDWR
      MODULE PROCEDURE FORCING_FIELDS_GET_HOST_DATA_RDWR
   END INTERFACE GET_HOST_DATA_RDWR
   INTERFACE GET_HOST_DATA_RDONLY
      MODULE PROCEDURE FORCING_FIELDS_GET_HOST_DATA_RDONLY
   END INTERFACE GET_HOST_DATA_RDONLY
   INTERFACE DELETE_DEVICE_DATA
      MODULE PROCEDURE FORCING_FIELDS_DELETE_DEVICE_DATA
   END INTERFACE DELETE_DEVICE_DATA

   TYPE FORCING_FIELDS
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_UWND => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UWND(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_UWND => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UWND(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VWND => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VWND(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VWND => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VWND(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_AIRD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: AIRD(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_AIRD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: AIRD(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSTAR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSTAR(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSTAR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSTAR(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CICOVER => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CICOVER(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CICOVER => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CICOVER(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CITHICK => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CITHICK(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CITHICK => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CITHICK(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_LKFR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: LKFR(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_LKFR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: LKFR(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_USTRA => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: USTRA(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_USTRA => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: USTRA(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VSTRA => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VSTRA(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VSTRA => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VSTRA(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_UCUR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UCUR(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_UCUR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UCUR(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VCUR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VCUR(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VCUR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VCUR(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSWAVE => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSWAVE(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSWAVE => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSWAVE(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WDWAVE => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WDWAVE(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WDWAVE => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WDWAVE(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_UFRIC => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UFRIC(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_UFRIC => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: UFRIC(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUW => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUW(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUW => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUW(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUWDIR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUWDIR(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUWDIR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUWDIR(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_Z0M => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0M(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_Z0M => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0M(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_Z0B => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0B(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_Z0B => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: Z0B(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CHRNCK => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CHRNCK(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CHRNCK => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CHRNCK(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_XLON => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XLON(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_XLON => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XLON(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_YLAT => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: YLAT(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_YLAT => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: YLAT(:,:) => NULL()
#endif
      LOGICAL :: LALLOC = .FALSE.
   END TYPE FORCING_FIELDS

   PUBLIC :: FORCING_FIELDS, ALLOC, DEALLOC
   PUBLIC :: SYNC_DEVICE_RDWR, SYNC_DEVICE_RDONLY, SYNC_HOST_RDWR, SYNC_HOST_RDONLY, GET_DEVICE_DATA_RDWR, GET_DEVICE_DATA_RDONLY,&
       & GET_DEVICE_DATA_WRONLY, GET_HOST_DATA_RDWR, GET_HOST_DATA_RDONLY, DELETE_DEVICE_DATA

   CONTAINS

   SUBROUTINE FORCING_FIELDS_ALLOC(SELF, UBOUNDS, LBOUNDS)
      TYPE(FORCING_FIELDS) :: SELF
      INTEGER(KIND=JWIM), INTENT(IN) :: UBOUNDS(2)
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: LBOUNDS(2)
      INTEGER(KIND=JWIM) :: LLBOUNDS(2)

      LLBOUNDS(:) = 1
      IF(PRESENT(LBOUNDS)) LLBOUNDS = LBOUNDS

      CALL FIELD_NEW(SELF%F_UWND, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_UWND%GET_HOST_DATA_RDWR(SELF%UWND)
      SELF%UWND(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_VWND, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_VWND%GET_HOST_DATA_RDWR(SELF%VWND)
      SELF%VWND(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_AIRD, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_AIRD%GET_HOST_DATA_RDWR(SELF%AIRD)
      SELF%AIRD(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_WSTAR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_WSTAR%GET_HOST_DATA_RDWR(SELF%WSTAR)
      SELF%WSTAR(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_CICOVER, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_CICOVER%GET_HOST_DATA_RDWR(SELF%CICOVER)
      SELF%CICOVER(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_CITHICK, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_CITHICK%GET_HOST_DATA_RDWR(SELF%CITHICK)
      SELF%CITHICK(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_LKFR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_LKFR%GET_HOST_DATA_RDWR(SELF%LKFR)
      SELF%LKFR(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_USTRA, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_USTRA%GET_HOST_DATA_RDWR(SELF%USTRA)
      SELF%USTRA(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_VSTRA, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_VSTRA%GET_HOST_DATA_RDWR(SELF%VSTRA)
      SELF%VSTRA(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_UCUR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_UCUR%GET_HOST_DATA_RDWR(SELF%UCUR)
      SELF%UCUR(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_VCUR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_VCUR%GET_HOST_DATA_RDWR(SELF%VCUR)
      SELF%VCUR(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_WSWAVE, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_WSWAVE%GET_HOST_DATA_RDWR(SELF%WSWAVE)
      SELF%WSWAVE(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_WDWAVE, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_WDWAVE%GET_HOST_DATA_RDWR(SELF%WDWAVE)
      SELF%WDWAVE(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_UFRIC, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_UFRIC%GET_HOST_DATA_RDWR(SELF%UFRIC)
      SELF%UFRIC(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUW, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUW%GET_HOST_DATA_RDWR(SELF%TAUW)
      SELF%TAUW(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUWDIR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUWDIR%GET_HOST_DATA_RDWR(SELF%TAUWDIR)
      SELF%TAUWDIR(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_Z0M, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_Z0M%GET_HOST_DATA_RDWR(SELF%Z0M)
      SELF%Z0M(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_Z0B, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_Z0B%GET_HOST_DATA_RDWR(SELF%Z0B)
      SELF%Z0B(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_CHRNCK, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_CHRNCK%GET_HOST_DATA_RDWR(SELF%CHRNCK)
      SELF%CHRNCK(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_XLON, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_XLON%GET_HOST_DATA_RDWR(SELF%XLON)
      SELF%XLON(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_YLAT, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_YLAT%GET_HOST_DATA_RDWR(SELF%YLAT)
      SELF%YLAT(:,:) = 0_JWRB

      SELF%LALLOC = .TRUE.
#ifdef OMPGPU
      !$omp target enter data map(to:SELF)
#else
      !$acc enter data copyin(SELF)
#endif

   END SUBROUTINE FORCING_FIELDS_ALLOC

   SUBROUTINE FORCING_FIELDS_DEALLOC(SELF)
      TYPE(FORCING_FIELDS) :: SELF

      NULLIFY(SELF%UWND)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%UWND)
#else
      !$acc exit data detach(SELF%UWND)
#endif
      CALL FIELD_DELETE(SELF%F_UWND)
      NULLIFY(SELF%F_UWND)
      NULLIFY(SELF%VWND)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%VWND)
#else
      !$acc exit data detach(SELF%VWND)
#endif
      CALL FIELD_DELETE(SELF%F_VWND)
      NULLIFY(SELF%F_VWND)
      NULLIFY(SELF%AIRD)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%AIRD)
#else
      !$acc exit data detach(SELF%AIRD)
#endif
      CALL FIELD_DELETE(SELF%F_AIRD)
      NULLIFY(SELF%F_AIRD)
      NULLIFY(SELF%WSTAR)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%WSTAR)
#else
      !$acc exit data detach(SELF%WSTAR)
#endif
      CALL FIELD_DELETE(SELF%F_WSTAR)
      NULLIFY(SELF%F_WSTAR)
      NULLIFY(SELF%CICOVER)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%CICOVER)
#else
      !$acc exit data detach(SELF%CICOVER)
#endif
      CALL FIELD_DELETE(SELF%F_CICOVER)
      NULLIFY(SELF%F_CICOVER)
      NULLIFY(SELF%CITHICK)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%CITHICK)
#else
      !$acc exit data detach(SELF%CITHICK)
#endif
      CALL FIELD_DELETE(SELF%F_CITHICK)
      NULLIFY(SELF%F_CITHICK)
      NULLIFY(SELF%LKFR)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%LKFR)
#else
      !$acc exit data detach(SELF%LKFR)
#endif
      CALL FIELD_DELETE(SELF%F_LKFR)
      NULLIFY(SELF%F_LKFR)
      NULLIFY(SELF%USTRA)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%USTRA)
#else
      !$acc exit data detach(SELF%USTRA)
#endif
      CALL FIELD_DELETE(SELF%F_USTRA)
      NULLIFY(SELF%F_USTRA)
      NULLIFY(SELF%VSTRA)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%VSTRA)
#else
      !$acc exit data detach(SELF%VSTRA)
#endif
      CALL FIELD_DELETE(SELF%F_VSTRA)
      NULLIFY(SELF%F_VSTRA)
      NULLIFY(SELF%UCUR)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%UCUR)
#else
      !$acc exit data detach(SELF%UCUR)
#endif
      CALL FIELD_DELETE(SELF%F_UCUR)
      NULLIFY(SELF%F_UCUR)
      NULLIFY(SELF%VCUR)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%VCUR)
#else
      !$acc exit data detach(SELF%VCUR)
#endif
      CALL FIELD_DELETE(SELF%F_VCUR)
      NULLIFY(SELF%F_VCUR)
      NULLIFY(SELF%WSWAVE)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%WSWAVE)
#else
      !$acc exit data detach(SELF%WSWAVE)
#endif
      CALL FIELD_DELETE(SELF%F_WSWAVE)
      NULLIFY(SELF%F_WSWAVE)
      NULLIFY(SELF%WDWAVE)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%WDWAVE)
#else
      !$acc exit data detach(SELF%WDWAVE)
#endif
      CALL FIELD_DELETE(SELF%F_WDWAVE)
      NULLIFY(SELF%F_WDWAVE)
      NULLIFY(SELF%UFRIC)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%UFRIC)
#else
      !$acc exit data detach(SELF%UFRIC)
#endif
      CALL FIELD_DELETE(SELF%F_UFRIC)
      NULLIFY(SELF%F_UFRIC)
      NULLIFY(SELF%TAUW)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUW)
#else
      !$acc exit data detach(SELF%TAUW)
#endif
      CALL FIELD_DELETE(SELF%F_TAUW)
      NULLIFY(SELF%F_TAUW)
      NULLIFY(SELF%TAUWDIR)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUWDIR)
#else
      !$acc exit data detach(SELF%TAUWDIR)
#endif
      CALL FIELD_DELETE(SELF%F_TAUWDIR)
      NULLIFY(SELF%F_TAUWDIR)
      NULLIFY(SELF%Z0M)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%Z0M)
#else
      !$acc exit data detach(SELF%Z0M)
#endif
      CALL FIELD_DELETE(SELF%F_Z0M)
      NULLIFY(SELF%F_Z0M)
      NULLIFY(SELF%Z0B)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%Z0B)
#else
      !$acc exit data detach(SELF%Z0B)
#endif
      CALL FIELD_DELETE(SELF%F_Z0B)
      NULLIFY(SELF%F_Z0B)
      NULLIFY(SELF%CHRNCK)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%CHRNCK)
#else
      !$acc exit data detach(SELF%CHRNCK)
#endif
      CALL FIELD_DELETE(SELF%F_CHRNCK)
      NULLIFY(SELF%F_CHRNCK)
      NULLIFY(SELF%XLON)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%XLON)
#else
      !$acc exit data detach(SELF%XLON)
#endif
      CALL FIELD_DELETE(SELF%F_XLON)
      NULLIFY(SELF%F_XLON)
      NULLIFY(SELF%YLAT)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%YLAT)
#else
      !$acc exit data detach(SELF%YLAT)
#endif
      CALL FIELD_DELETE(SELF%F_YLAT)
      NULLIFY(SELF%F_YLAT)

      SELF%LALLOC = .FALSE.
#ifdef OMPGPU
      !$omp target exit data map(delete:SELF)
#else
      !$acc exit data delete(SELF)
#endif

   END SUBROUTINE FORCING_FIELDS_DEALLOC

   SUBROUTINE FORCING_FIELDS_GET_DEVICE_DATA_RDWR(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR, VCUR,&
       & WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
             CALL SELF%F_UWND%GET_DEVICE_DATA_RDWR(SELF%UWND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UWND)
#else
             !$acc enter data attach(SELF%UWND)
#endif
         ENDIF
         IF(L_VWND)THEN
             CALL SELF%F_VWND%GET_DEVICE_DATA_RDWR(SELF%VWND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VWND)
#else
             !$acc enter data attach(SELF%VWND)
#endif
         ENDIF
         IF(L_AIRD)THEN
             CALL SELF%F_AIRD%GET_DEVICE_DATA_RDWR(SELF%AIRD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%AIRD)
#else
             !$acc enter data attach(SELF%AIRD)
#endif
         ENDIF
         IF(L_WSTAR)THEN
             CALL SELF%F_WSTAR%GET_DEVICE_DATA_RDWR(SELF%WSTAR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSTAR)
#else
             !$acc enter data attach(SELF%WSTAR)
#endif
         ENDIF
         IF(L_CICOVER)THEN
             CALL SELF%F_CICOVER%GET_DEVICE_DATA_RDWR(SELF%CICOVER)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CICOVER)
#else
             !$acc enter data attach(SELF%CICOVER)
#endif
         ENDIF
         IF(L_CITHICK)THEN
             CALL SELF%F_CITHICK%GET_DEVICE_DATA_RDWR(SELF%CITHICK)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CITHICK)
#else
             !$acc enter data attach(SELF%CITHICK)
#endif
         ENDIF
         IF(L_LKFR)THEN
             CALL SELF%F_LKFR%GET_DEVICE_DATA_RDWR(SELF%LKFR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%LKFR)
#else
             !$acc enter data attach(SELF%LKFR)
#endif
         ENDIF
         IF(L_USTRA)THEN
             CALL SELF%F_USTRA%GET_DEVICE_DATA_RDWR(SELF%USTRA)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%USTRA)
#else
             !$acc enter data attach(SELF%USTRA)
#endif
         ENDIF
         IF(L_VSTRA)THEN
             CALL SELF%F_VSTRA%GET_DEVICE_DATA_RDWR(SELF%VSTRA)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VSTRA)
#else
             !$acc enter data attach(SELF%VSTRA)
#endif
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%GET_DEVICE_DATA_RDWR(SELF%UCUR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UCUR)
#else
             !$acc enter data attach(SELF%UCUR)
#endif
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%GET_DEVICE_DATA_RDWR(SELF%VCUR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VCUR)
#else
             !$acc enter data attach(SELF%VCUR)
#endif
         ENDIF
         IF(L_WSWAVE)THEN
             CALL SELF%F_WSWAVE%GET_DEVICE_DATA_RDWR(SELF%WSWAVE)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSWAVE)
#else
             !$acc enter data attach(SELF%WSWAVE)
#endif
         ENDIF
         IF(L_WDWAVE)THEN
             CALL SELF%F_WDWAVE%GET_DEVICE_DATA_RDWR(SELF%WDWAVE)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WDWAVE)
#else
             !$acc enter data attach(SELF%WDWAVE)
#endif
         ENDIF
         IF(L_UFRIC)THEN
             CALL SELF%F_UFRIC%GET_DEVICE_DATA_RDWR(SELF%UFRIC)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UFRIC)
#else
             !$acc enter data attach(SELF%UFRIC)
#endif
         ENDIF
         IF(L_TAUW)THEN
             CALL SELF%F_TAUW%GET_DEVICE_DATA_RDWR(SELF%TAUW)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUW)
#else
             !$acc enter data attach(SELF%TAUW)
#endif
         ENDIF
         IF(L_TAUWDIR)THEN
             CALL SELF%F_TAUWDIR%GET_DEVICE_DATA_RDWR(SELF%TAUWDIR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUWDIR)
#else
             !$acc enter data attach(SELF%TAUWDIR)
#endif
         ENDIF
         IF(L_Z0M)THEN
             CALL SELF%F_Z0M%GET_DEVICE_DATA_RDWR(SELF%Z0M)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%Z0M)
#else
             !$acc enter data attach(SELF%Z0M)
#endif
         ENDIF
         IF(L_Z0B)THEN
             CALL SELF%F_Z0B%GET_DEVICE_DATA_RDWR(SELF%Z0B)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%Z0B)
#else
             !$acc enter data attach(SELF%Z0B)
#endif
         ENDIF
         IF(L_CHRNCK)THEN
             CALL SELF%F_CHRNCK%GET_DEVICE_DATA_RDWR(SELF%CHRNCK)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CHRNCK)
#else
             !$acc enter data attach(SELF%CHRNCK)
#endif
         ENDIF
         IF(L_XLON)THEN
             CALL SELF%F_XLON%GET_DEVICE_DATA_RDWR(SELF%XLON)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%XLON)
#else
             !$acc enter data attach(SELF%XLON)
#endif
         ENDIF
         IF(L_YLAT)THEN
             CALL SELF%F_YLAT%GET_DEVICE_DATA_RDWR(SELF%YLAT)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%YLAT)
#else
             !$acc enter data attach(SELF%YLAT)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_UWND%GET_DEVICE_DATA_RDWR(SELF%UWND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UWND)
#else
         !$acc enter data attach(SELF%UWND)
#endif
         CALL SELF%F_VWND%GET_DEVICE_DATA_RDWR(SELF%VWND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VWND)
#else
         !$acc enter data attach(SELF%VWND)
#endif
         CALL SELF%F_AIRD%GET_DEVICE_DATA_RDWR(SELF%AIRD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%AIRD)
#else
         !$acc enter data attach(SELF%AIRD)
#endif
         CALL SELF%F_WSTAR%GET_DEVICE_DATA_RDWR(SELF%WSTAR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSTAR)
#else
         !$acc enter data attach(SELF%WSTAR)
#endif
         CALL SELF%F_CICOVER%GET_DEVICE_DATA_RDWR(SELF%CICOVER)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CICOVER)
#else
         !$acc enter data attach(SELF%CICOVER)
#endif
         CALL SELF%F_CITHICK%GET_DEVICE_DATA_RDWR(SELF%CITHICK)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CITHICK)
#else
         !$acc enter data attach(SELF%CITHICK)
#endif
         CALL SELF%F_LKFR%GET_DEVICE_DATA_RDWR(SELF%LKFR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%LKFR)
#else
         !$acc enter data attach(SELF%LKFR)
#endif
         CALL SELF%F_USTRA%GET_DEVICE_DATA_RDWR(SELF%USTRA)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%USTRA)
#else
         !$acc enter data attach(SELF%USTRA)
#endif
         CALL SELF%F_VSTRA%GET_DEVICE_DATA_RDWR(SELF%VSTRA)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VSTRA)
#else
         !$acc enter data attach(SELF%VSTRA)
#endif
         CALL SELF%F_UCUR%GET_DEVICE_DATA_RDWR(SELF%UCUR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UCUR)
#else
         !$acc enter data attach(SELF%UCUR)
#endif
         CALL SELF%F_VCUR%GET_DEVICE_DATA_RDWR(SELF%VCUR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VCUR)
#else
         !$acc enter data attach(SELF%VCUR)
#endif
         CALL SELF%F_WSWAVE%GET_DEVICE_DATA_RDWR(SELF%WSWAVE)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSWAVE)
#else
         !$acc enter data attach(SELF%WSWAVE)
#endif
         CALL SELF%F_WDWAVE%GET_DEVICE_DATA_RDWR(SELF%WDWAVE)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WDWAVE)
#else
         !$acc enter data attach(SELF%WDWAVE)
#endif
         CALL SELF%F_UFRIC%GET_DEVICE_DATA_RDWR(SELF%UFRIC)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UFRIC)
#else
         !$acc enter data attach(SELF%UFRIC)
#endif
         CALL SELF%F_TAUW%GET_DEVICE_DATA_RDWR(SELF%TAUW)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUW)
#else
         !$acc enter data attach(SELF%TAUW)
#endif
         CALL SELF%F_TAUWDIR%GET_DEVICE_DATA_RDWR(SELF%TAUWDIR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUWDIR)
#else
         !$acc enter data attach(SELF%TAUWDIR)
#endif
         CALL SELF%F_Z0M%GET_DEVICE_DATA_RDWR(SELF%Z0M)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%Z0M)
#else
         !$acc enter data attach(SELF%Z0M)
#endif
         CALL SELF%F_Z0B%GET_DEVICE_DATA_RDWR(SELF%Z0B)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%Z0B)
#else
         !$acc enter data attach(SELF%Z0B)
#endif
         CALL SELF%F_CHRNCK%GET_DEVICE_DATA_RDWR(SELF%CHRNCK)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CHRNCK)
#else
         !$acc enter data attach(SELF%CHRNCK)
#endif
         CALL SELF%F_XLON%GET_DEVICE_DATA_RDWR(SELF%XLON)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%XLON)
#else
         !$acc enter data attach(SELF%XLON)
#endif
         CALL SELF%F_YLAT%GET_DEVICE_DATA_RDWR(SELF%YLAT)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%YLAT)
#else
         !$acc enter data attach(SELF%YLAT)
#endif
      ENDIF

   END SUBROUTINE FORCING_FIELDS_GET_DEVICE_DATA_RDWR

   SUBROUTINE FORCING_FIELDS_GET_DEVICE_DATA_WRONLY(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR,&
       & VCUR, WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
             CALL SELF%F_UWND%GET_DEVICE_DATA_WRONLY(SELF%UWND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UWND)
#else
             !$acc enter data attach(SELF%UWND)
#endif
         ENDIF
         IF(L_VWND)THEN
             CALL SELF%F_VWND%GET_DEVICE_DATA_WRONLY(SELF%VWND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VWND)
#else
             !$acc enter data attach(SELF%VWND)
#endif
         ENDIF
         IF(L_AIRD)THEN
             CALL SELF%F_AIRD%GET_DEVICE_DATA_WRONLY(SELF%AIRD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%AIRD)
#else
             !$acc enter data attach(SELF%AIRD)
#endif
         ENDIF
         IF(L_WSTAR)THEN
             CALL SELF%F_WSTAR%GET_DEVICE_DATA_WRONLY(SELF%WSTAR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSTAR)
#else
             !$acc enter data attach(SELF%WSTAR)
#endif
         ENDIF
         IF(L_CICOVER)THEN
             CALL SELF%F_CICOVER%GET_DEVICE_DATA_WRONLY(SELF%CICOVER)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CICOVER)
#else
             !$acc enter data attach(SELF%CICOVER)
#endif
         ENDIF
         IF(L_CITHICK)THEN
             CALL SELF%F_CITHICK%GET_DEVICE_DATA_WRONLY(SELF%CITHICK)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CITHICK)
#else
             !$acc enter data attach(SELF%CITHICK)
#endif
         ENDIF
         IF(L_LKFR)THEN
             CALL SELF%F_LKFR%GET_DEVICE_DATA_WRONLY(SELF%LKFR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%LKFR)
#else
             !$acc enter data attach(SELF%LKFR)
#endif
         ENDIF
         IF(L_USTRA)THEN
             CALL SELF%F_USTRA%GET_DEVICE_DATA_WRONLY(SELF%USTRA)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%USTRA)
#else
             !$acc enter data attach(SELF%USTRA)
#endif
         ENDIF
         IF(L_VSTRA)THEN
             CALL SELF%F_VSTRA%GET_DEVICE_DATA_WRONLY(SELF%VSTRA)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VSTRA)
#else
             !$acc enter data attach(SELF%VSTRA)
#endif
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%GET_DEVICE_DATA_WRONLY(SELF%UCUR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UCUR)
#else
             !$acc enter data attach(SELF%UCUR)
#endif
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%GET_DEVICE_DATA_WRONLY(SELF%VCUR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VCUR)
#else
             !$acc enter data attach(SELF%VCUR)
#endif
         ENDIF
         IF(L_WSWAVE)THEN
             CALL SELF%F_WSWAVE%GET_DEVICE_DATA_WRONLY(SELF%WSWAVE)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSWAVE)
#else
             !$acc enter data attach(SELF%WSWAVE)
#endif
         ENDIF
         IF(L_WDWAVE)THEN
             CALL SELF%F_WDWAVE%GET_DEVICE_DATA_WRONLY(SELF%WDWAVE)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WDWAVE)
#else
             !$acc enter data attach(SELF%WDWAVE)
#endif
         ENDIF
         IF(L_UFRIC)THEN
             CALL SELF%F_UFRIC%GET_DEVICE_DATA_WRONLY(SELF%UFRIC)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UFRIC)
#else
             !$acc enter data attach(SELF%UFRIC)
#endif
         ENDIF
         IF(L_TAUW)THEN
             CALL SELF%F_TAUW%GET_DEVICE_DATA_WRONLY(SELF%TAUW)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUW)
#else
             !$acc enter data attach(SELF%TAUW)
#endif
         ENDIF
         IF(L_TAUWDIR)THEN
             CALL SELF%F_TAUWDIR%GET_DEVICE_DATA_WRONLY(SELF%TAUWDIR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUWDIR)
#else
             !$acc enter data attach(SELF%TAUWDIR)
#endif
         ENDIF
         IF(L_Z0M)THEN
             CALL SELF%F_Z0M%GET_DEVICE_DATA_WRONLY(SELF%Z0M)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%Z0M)
#else
             !$acc enter data attach(SELF%Z0M)
#endif
         ENDIF
         IF(L_Z0B)THEN
             CALL SELF%F_Z0B%GET_DEVICE_DATA_WRONLY(SELF%Z0B)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%Z0B)
#else
             !$acc enter data attach(SELF%Z0B)
#endif
         ENDIF
         IF(L_CHRNCK)THEN
             CALL SELF%F_CHRNCK%GET_DEVICE_DATA_WRONLY(SELF%CHRNCK)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CHRNCK)
#else
             !$acc enter data attach(SELF%CHRNCK)
#endif
         ENDIF
         IF(L_XLON)THEN
             CALL SELF%F_XLON%GET_DEVICE_DATA_WRONLY(SELF%XLON)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%XLON)
#else
             !$acc enter data attach(SELF%XLON)
#endif
         ENDIF
         IF(L_YLAT)THEN
             CALL SELF%F_YLAT%GET_DEVICE_DATA_WRONLY(SELF%YLAT)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%YLAT)
#else
             !$acc enter data attach(SELF%YLAT)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_UWND%GET_DEVICE_DATA_WRONLY(SELF%UWND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UWND)
#else
         !$acc enter data attach(SELF%UWND)
#endif
         CALL SELF%F_VWND%GET_DEVICE_DATA_WRONLY(SELF%VWND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VWND)
#else
         !$acc enter data attach(SELF%VWND)
#endif
         CALL SELF%F_AIRD%GET_DEVICE_DATA_WRONLY(SELF%AIRD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%AIRD)
#else
         !$acc enter data attach(SELF%AIRD)
#endif
         CALL SELF%F_WSTAR%GET_DEVICE_DATA_WRONLY(SELF%WSTAR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSTAR)
#else
         !$acc enter data attach(SELF%WSTAR)
#endif
         CALL SELF%F_CICOVER%GET_DEVICE_DATA_WRONLY(SELF%CICOVER)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CICOVER)
#else
         !$acc enter data attach(SELF%CICOVER)
#endif
         CALL SELF%F_CITHICK%GET_DEVICE_DATA_WRONLY(SELF%CITHICK)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CITHICK)
#else
         !$acc enter data attach(SELF%CITHICK)
#endif
         CALL SELF%F_LKFR%GET_DEVICE_DATA_WRONLY(SELF%LKFR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%LKFR)
#else
         !$acc enter data attach(SELF%LKFR)
#endif
         CALL SELF%F_USTRA%GET_DEVICE_DATA_WRONLY(SELF%USTRA)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%USTRA)
#else
         !$acc enter data attach(SELF%USTRA)
#endif
         CALL SELF%F_VSTRA%GET_DEVICE_DATA_WRONLY(SELF%VSTRA)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VSTRA)
#else
         !$acc enter data attach(SELF%VSTRA)
#endif
         CALL SELF%F_UCUR%GET_DEVICE_DATA_WRONLY(SELF%UCUR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UCUR)
#else
         !$acc enter data attach(SELF%UCUR)
#endif
         CALL SELF%F_VCUR%GET_DEVICE_DATA_WRONLY(SELF%VCUR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VCUR)
#else
         !$acc enter data attach(SELF%VCUR)
#endif
         CALL SELF%F_WSWAVE%GET_DEVICE_DATA_WRONLY(SELF%WSWAVE)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSWAVE)
#else
         !$acc enter data attach(SELF%WSWAVE)
#endif
         CALL SELF%F_WDWAVE%GET_DEVICE_DATA_WRONLY(SELF%WDWAVE)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WDWAVE)
#else
         !$acc enter data attach(SELF%WDWAVE)
#endif
         CALL SELF%F_UFRIC%GET_DEVICE_DATA_WRONLY(SELF%UFRIC)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UFRIC)
#else
         !$acc enter data attach(SELF%UFRIC)
#endif
         CALL SELF%F_TAUW%GET_DEVICE_DATA_WRONLY(SELF%TAUW)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUW)
#else
         !$acc enter data attach(SELF%TAUW)
#endif
         CALL SELF%F_TAUWDIR%GET_DEVICE_DATA_WRONLY(SELF%TAUWDIR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUWDIR)
#else
         !$acc enter data attach(SELF%TAUWDIR)
#endif
         CALL SELF%F_Z0M%GET_DEVICE_DATA_WRONLY(SELF%Z0M)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%Z0M)
#else
         !$acc enter data attach(SELF%Z0M)
#endif
         CALL SELF%F_Z0B%GET_DEVICE_DATA_WRONLY(SELF%Z0B)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%Z0B)
#else
         !$acc enter data attach(SELF%Z0B)
#endif
         CALL SELF%F_CHRNCK%GET_DEVICE_DATA_WRONLY(SELF%CHRNCK)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CHRNCK)
#else
         !$acc enter data attach(SELF%CHRNCK)
#endif
         CALL SELF%F_XLON%GET_DEVICE_DATA_WRONLY(SELF%XLON)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%XLON)
#else
         !$acc enter data attach(SELF%XLON)
#endif
         CALL SELF%F_YLAT%GET_DEVICE_DATA_WRONLY(SELF%YLAT)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%YLAT)
#else
         !$acc enter data attach(SELF%YLAT)
#endif
      ENDIF

   END SUBROUTINE FORCING_FIELDS_GET_DEVICE_DATA_WRONLY

   SUBROUTINE FORCING_FIELDS_GET_DEVICE_DATA_RDONLY(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR,&
       & VCUR, WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
             CALL SELF%F_UWND%GET_DEVICE_DATA_RDONLY(SELF%UWND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UWND)
#else
             !$acc enter data attach(SELF%UWND)
#endif
         ENDIF
         IF(L_VWND)THEN
             CALL SELF%F_VWND%GET_DEVICE_DATA_RDONLY(SELF%VWND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VWND)
#else
             !$acc enter data attach(SELF%VWND)
#endif
         ENDIF
         IF(L_AIRD)THEN
             CALL SELF%F_AIRD%GET_DEVICE_DATA_RDONLY(SELF%AIRD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%AIRD)
#else
             !$acc enter data attach(SELF%AIRD)
#endif
         ENDIF
         IF(L_WSTAR)THEN
             CALL SELF%F_WSTAR%GET_DEVICE_DATA_RDONLY(SELF%WSTAR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSTAR)
#else
             !$acc enter data attach(SELF%WSTAR)
#endif
         ENDIF
         IF(L_CICOVER)THEN
             CALL SELF%F_CICOVER%GET_DEVICE_DATA_RDONLY(SELF%CICOVER)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CICOVER)
#else
             !$acc enter data attach(SELF%CICOVER)
#endif
         ENDIF
         IF(L_CITHICK)THEN
             CALL SELF%F_CITHICK%GET_DEVICE_DATA_RDONLY(SELF%CITHICK)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CITHICK)
#else
             !$acc enter data attach(SELF%CITHICK)
#endif
         ENDIF
         IF(L_LKFR)THEN
             CALL SELF%F_LKFR%GET_DEVICE_DATA_RDONLY(SELF%LKFR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%LKFR)
#else
             !$acc enter data attach(SELF%LKFR)
#endif
         ENDIF
         IF(L_USTRA)THEN
             CALL SELF%F_USTRA%GET_DEVICE_DATA_RDONLY(SELF%USTRA)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%USTRA)
#else
             !$acc enter data attach(SELF%USTRA)
#endif
         ENDIF
         IF(L_VSTRA)THEN
             CALL SELF%F_VSTRA%GET_DEVICE_DATA_RDONLY(SELF%VSTRA)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VSTRA)
#else
             !$acc enter data attach(SELF%VSTRA)
#endif
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%GET_DEVICE_DATA_RDONLY(SELF%UCUR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UCUR)
#else
             !$acc enter data attach(SELF%UCUR)
#endif
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%GET_DEVICE_DATA_RDONLY(SELF%VCUR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VCUR)
#else
             !$acc enter data attach(SELF%VCUR)
#endif
         ENDIF
         IF(L_WSWAVE)THEN
             CALL SELF%F_WSWAVE%GET_DEVICE_DATA_RDONLY(SELF%WSWAVE)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSWAVE)
#else
             !$acc enter data attach(SELF%WSWAVE)
#endif
         ENDIF
         IF(L_WDWAVE)THEN
             CALL SELF%F_WDWAVE%GET_DEVICE_DATA_RDONLY(SELF%WDWAVE)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WDWAVE)
#else
             !$acc enter data attach(SELF%WDWAVE)
#endif
         ENDIF
         IF(L_UFRIC)THEN
             CALL SELF%F_UFRIC%GET_DEVICE_DATA_RDONLY(SELF%UFRIC)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%UFRIC)
#else
             !$acc enter data attach(SELF%UFRIC)
#endif
         ENDIF
         IF(L_TAUW)THEN
             CALL SELF%F_TAUW%GET_DEVICE_DATA_RDONLY(SELF%TAUW)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUW)
#else
             !$acc enter data attach(SELF%TAUW)
#endif
         ENDIF
         IF(L_TAUWDIR)THEN
             CALL SELF%F_TAUWDIR%GET_DEVICE_DATA_RDONLY(SELF%TAUWDIR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUWDIR)
#else
             !$acc enter data attach(SELF%TAUWDIR)
#endif
         ENDIF
         IF(L_Z0M)THEN
             CALL SELF%F_Z0M%GET_DEVICE_DATA_RDONLY(SELF%Z0M)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%Z0M)
#else
             !$acc enter data attach(SELF%Z0M)
#endif
         ENDIF
         IF(L_Z0B)THEN
             CALL SELF%F_Z0B%GET_DEVICE_DATA_RDONLY(SELF%Z0B)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%Z0B)
#else
             !$acc enter data attach(SELF%Z0B)
#endif
         ENDIF
         IF(L_CHRNCK)THEN
             CALL SELF%F_CHRNCK%GET_DEVICE_DATA_RDONLY(SELF%CHRNCK)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CHRNCK)
#else
             !$acc enter data attach(SELF%CHRNCK)
#endif
         ENDIF
         IF(L_XLON)THEN
             CALL SELF%F_XLON%GET_DEVICE_DATA_RDONLY(SELF%XLON)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%XLON)
#else
             !$acc enter data attach(SELF%XLON)
#endif
         ENDIF
         IF(L_YLAT)THEN
             CALL SELF%F_YLAT%GET_DEVICE_DATA_RDONLY(SELF%YLAT)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%YLAT)
#else
             !$acc enter data attach(SELF%YLAT)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_UWND%GET_DEVICE_DATA_RDONLY(SELF%UWND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UWND)
#else
         !$acc enter data attach(SELF%UWND)
#endif
         CALL SELF%F_VWND%GET_DEVICE_DATA_RDONLY(SELF%VWND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VWND)
#else
         !$acc enter data attach(SELF%VWND)
#endif
         CALL SELF%F_AIRD%GET_DEVICE_DATA_RDONLY(SELF%AIRD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%AIRD)
#else
         !$acc enter data attach(SELF%AIRD)
#endif
         CALL SELF%F_WSTAR%GET_DEVICE_DATA_RDONLY(SELF%WSTAR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSTAR)
#else
         !$acc enter data attach(SELF%WSTAR)
#endif
         CALL SELF%F_CICOVER%GET_DEVICE_DATA_RDONLY(SELF%CICOVER)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CICOVER)
#else
         !$acc enter data attach(SELF%CICOVER)
#endif
         CALL SELF%F_CITHICK%GET_DEVICE_DATA_RDONLY(SELF%CITHICK)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CITHICK)
#else
         !$acc enter data attach(SELF%CITHICK)
#endif
         CALL SELF%F_LKFR%GET_DEVICE_DATA_RDONLY(SELF%LKFR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%LKFR)
#else
         !$acc enter data attach(SELF%LKFR)
#endif
         CALL SELF%F_USTRA%GET_DEVICE_DATA_RDONLY(SELF%USTRA)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%USTRA)
#else
         !$acc enter data attach(SELF%USTRA)
#endif
         CALL SELF%F_VSTRA%GET_DEVICE_DATA_RDONLY(SELF%VSTRA)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VSTRA)
#else
         !$acc enter data attach(SELF%VSTRA)
#endif
         CALL SELF%F_UCUR%GET_DEVICE_DATA_RDONLY(SELF%UCUR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UCUR)
#else
         !$acc enter data attach(SELF%UCUR)
#endif
         CALL SELF%F_VCUR%GET_DEVICE_DATA_RDONLY(SELF%VCUR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VCUR)
#else
         !$acc enter data attach(SELF%VCUR)
#endif
         CALL SELF%F_WSWAVE%GET_DEVICE_DATA_RDONLY(SELF%WSWAVE)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSWAVE)
#else
         !$acc enter data attach(SELF%WSWAVE)
#endif
         CALL SELF%F_WDWAVE%GET_DEVICE_DATA_RDONLY(SELF%WDWAVE)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WDWAVE)
#else
         !$acc enter data attach(SELF%WDWAVE)
#endif
         CALL SELF%F_UFRIC%GET_DEVICE_DATA_RDONLY(SELF%UFRIC)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%UFRIC)
#else
         !$acc enter data attach(SELF%UFRIC)
#endif
         CALL SELF%F_TAUW%GET_DEVICE_DATA_RDONLY(SELF%TAUW)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUW)
#else
         !$acc enter data attach(SELF%TAUW)
#endif
         CALL SELF%F_TAUWDIR%GET_DEVICE_DATA_RDONLY(SELF%TAUWDIR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUWDIR)
#else
         !$acc enter data attach(SELF%TAUWDIR)
#endif
         CALL SELF%F_Z0M%GET_DEVICE_DATA_RDONLY(SELF%Z0M)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%Z0M)
#else
         !$acc enter data attach(SELF%Z0M)
#endif
         CALL SELF%F_Z0B%GET_DEVICE_DATA_RDONLY(SELF%Z0B)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%Z0B)
#else
         !$acc enter data attach(SELF%Z0B)
#endif
         CALL SELF%F_CHRNCK%GET_DEVICE_DATA_RDONLY(SELF%CHRNCK)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CHRNCK)
#else
         !$acc enter data attach(SELF%CHRNCK)
#endif
         CALL SELF%F_XLON%GET_DEVICE_DATA_RDONLY(SELF%XLON)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%XLON)
#else
         !$acc enter data attach(SELF%XLON)
#endif
         CALL SELF%F_YLAT%GET_DEVICE_DATA_RDONLY(SELF%YLAT)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%YLAT)
#else
         !$acc enter data attach(SELF%YLAT)
#endif
      ENDIF

   END SUBROUTINE FORCING_FIELDS_GET_DEVICE_DATA_RDONLY

   SUBROUTINE FORCING_FIELDS_SYNC_DEVICE_RDWR(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR, VCUR,&
       & WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT, QUEUE)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
             CALL SELF%F_UWND%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VWND)THEN
             CALL SELF%F_VWND%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_AIRD)THEN
             CALL SELF%F_AIRD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WSTAR)THEN
             CALL SELF%F_WSTAR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CICOVER)THEN
             CALL SELF%F_CICOVER%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CITHICK)THEN
             CALL SELF%F_CITHICK%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_LKFR)THEN
             CALL SELF%F_LKFR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_USTRA)THEN
             CALL SELF%F_USTRA%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTRA)THEN
             CALL SELF%F_VSTRA%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WSWAVE)THEN
             CALL SELF%F_WSWAVE%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WDWAVE)THEN
             CALL SELF%F_WDWAVE%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_UFRIC)THEN
             CALL SELF%F_UFRIC%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUW)THEN
             CALL SELF%F_TAUW%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUWDIR)THEN
             CALL SELF%F_TAUWDIR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0M)THEN
             CALL SELF%F_Z0M%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0B)THEN
             CALL SELF%F_Z0B%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CHRNCK)THEN
             CALL SELF%F_CHRNCK%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_XLON)THEN
             CALL SELF%F_XLON%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_YLAT)THEN
             CALL SELF%F_YLAT%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_UWND%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VWND%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_AIRD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WSTAR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CICOVER%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CITHICK%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_LKFR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_USTRA%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VSTRA%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WSWAVE%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WDWAVE%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_UFRIC%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUW%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUWDIR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_Z0M%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_Z0B%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CHRNCK%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_XLON%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_YLAT%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE FORCING_FIELDS_SYNC_DEVICE_RDWR

   SUBROUTINE FORCING_FIELDS_SYNC_DEVICE_RDONLY(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR, VCUR,&
       & WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT, QUEUE)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
             CALL SELF%F_UWND%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VWND)THEN
             CALL SELF%F_VWND%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_AIRD)THEN
             CALL SELF%F_AIRD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WSTAR)THEN
             CALL SELF%F_WSTAR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CICOVER)THEN
             CALL SELF%F_CICOVER%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CITHICK)THEN
             CALL SELF%F_CITHICK%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_LKFR)THEN
             CALL SELF%F_LKFR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_USTRA)THEN
             CALL SELF%F_USTRA%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTRA)THEN
             CALL SELF%F_VSTRA%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WSWAVE)THEN
             CALL SELF%F_WSWAVE%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WDWAVE)THEN
             CALL SELF%F_WDWAVE%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_UFRIC)THEN
             CALL SELF%F_UFRIC%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUW)THEN
             CALL SELF%F_TAUW%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUWDIR)THEN
             CALL SELF%F_TAUWDIR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0M)THEN
             CALL SELF%F_Z0M%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0B)THEN
             CALL SELF%F_Z0B%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CHRNCK)THEN
             CALL SELF%F_CHRNCK%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_XLON)THEN
             CALL SELF%F_XLON%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_YLAT)THEN
             CALL SELF%F_YLAT%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_UWND%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VWND%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_AIRD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WSTAR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CICOVER%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CITHICK%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_LKFR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_USTRA%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VSTRA%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WSWAVE%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WDWAVE%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_UFRIC%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUW%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUWDIR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_Z0M%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_Z0B%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CHRNCK%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_XLON%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_YLAT%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE FORCING_FIELDS_SYNC_DEVICE_RDONLY

   SUBROUTINE FORCING_FIELDS_GET_HOST_DATA_RDONLY(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR, VCUR,&
       & WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%UWND)
#else
             !$acc exit data detach(SELF%UWND)
#endif
             CALL SELF%F_UWND%GET_HOST_DATA_RDONLY(SELF%UWND)
         ENDIF
         IF(L_VWND)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VWND)
#else
             !$acc exit data detach(SELF%VWND)
#endif
             CALL SELF%F_VWND%GET_HOST_DATA_RDONLY(SELF%VWND)
         ENDIF
         IF(L_AIRD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%AIRD)
#else
             !$acc exit data detach(SELF%AIRD)
#endif
             CALL SELF%F_AIRD%GET_HOST_DATA_RDONLY(SELF%AIRD)
         ENDIF
         IF(L_WSTAR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSTAR)
#else
             !$acc exit data detach(SELF%WSTAR)
#endif
             CALL SELF%F_WSTAR%GET_HOST_DATA_RDONLY(SELF%WSTAR)
         ENDIF
         IF(L_CICOVER)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CICOVER)
#else
             !$acc exit data detach(SELF%CICOVER)
#endif
             CALL SELF%F_CICOVER%GET_HOST_DATA_RDONLY(SELF%CICOVER)
         ENDIF
         IF(L_CITHICK)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CITHICK)
#else
             !$acc exit data detach(SELF%CITHICK)
#endif
             CALL SELF%F_CITHICK%GET_HOST_DATA_RDONLY(SELF%CITHICK)
         ENDIF
         IF(L_LKFR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%LKFR)
#else
             !$acc exit data detach(SELF%LKFR)
#endif
             CALL SELF%F_LKFR%GET_HOST_DATA_RDONLY(SELF%LKFR)
         ENDIF
         IF(L_USTRA)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%USTRA)
#else
             !$acc exit data detach(SELF%USTRA)
#endif
             CALL SELF%F_USTRA%GET_HOST_DATA_RDONLY(SELF%USTRA)
         ENDIF
         IF(L_VSTRA)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VSTRA)
#else
             !$acc exit data detach(SELF%VSTRA)
#endif
             CALL SELF%F_VSTRA%GET_HOST_DATA_RDONLY(SELF%VSTRA)
         ENDIF
         IF(L_UCUR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%UCUR)
#else
             !$acc exit data detach(SELF%UCUR)
#endif
             CALL SELF%F_UCUR%GET_HOST_DATA_RDONLY(SELF%UCUR)
         ENDIF
         IF(L_VCUR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VCUR)
#else
             !$acc exit data detach(SELF%VCUR)
#endif
             CALL SELF%F_VCUR%GET_HOST_DATA_RDONLY(SELF%VCUR)
         ENDIF
         IF(L_WSWAVE)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSWAVE)
#else
             !$acc exit data detach(SELF%WSWAVE)
#endif
             CALL SELF%F_WSWAVE%GET_HOST_DATA_RDONLY(SELF%WSWAVE)
         ENDIF
         IF(L_WDWAVE)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WDWAVE)
#else
             !$acc exit data detach(SELF%WDWAVE)
#endif
             CALL SELF%F_WDWAVE%GET_HOST_DATA_RDONLY(SELF%WDWAVE)
         ENDIF
         IF(L_UFRIC)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%UFRIC)
#else
             !$acc exit data detach(SELF%UFRIC)
#endif
             CALL SELF%F_UFRIC%GET_HOST_DATA_RDONLY(SELF%UFRIC)
         ENDIF
         IF(L_TAUW)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUW)
#else
             !$acc exit data detach(SELF%TAUW)
#endif
             CALL SELF%F_TAUW%GET_HOST_DATA_RDONLY(SELF%TAUW)
         ENDIF
         IF(L_TAUWDIR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUWDIR)
#else
             !$acc exit data detach(SELF%TAUWDIR)
#endif
             CALL SELF%F_TAUWDIR%GET_HOST_DATA_RDONLY(SELF%TAUWDIR)
         ENDIF
         IF(L_Z0M)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%Z0M)
#else
             !$acc exit data detach(SELF%Z0M)
#endif
             CALL SELF%F_Z0M%GET_HOST_DATA_RDONLY(SELF%Z0M)
         ENDIF
         IF(L_Z0B)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%Z0B)
#else
             !$acc exit data detach(SELF%Z0B)
#endif
             CALL SELF%F_Z0B%GET_HOST_DATA_RDONLY(SELF%Z0B)
         ENDIF
         IF(L_CHRNCK)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CHRNCK)
#else
             !$acc exit data detach(SELF%CHRNCK)
#endif
             CALL SELF%F_CHRNCK%GET_HOST_DATA_RDONLY(SELF%CHRNCK)
         ENDIF
         IF(L_XLON)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%XLON)
#else
             !$acc exit data detach(SELF%XLON)
#endif
             CALL SELF%F_XLON%GET_HOST_DATA_RDONLY(SELF%XLON)
         ENDIF
         IF(L_YLAT)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%YLAT)
#else
             !$acc exit data detach(SELF%YLAT)
#endif
             CALL SELF%F_YLAT%GET_HOST_DATA_RDONLY(SELF%YLAT)
         ENDIF
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%UWND)
#else
         !$acc exit data detach(SELF%UWND)
#endif
         CALL SELF%F_UWND%GET_HOST_DATA_RDONLY(SELF%UWND)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VWND)
#else
         !$acc exit data detach(SELF%VWND)
#endif
         CALL SELF%F_VWND%GET_HOST_DATA_RDONLY(SELF%VWND)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%AIRD)
#else
         !$acc exit data detach(SELF%AIRD)
#endif
         CALL SELF%F_AIRD%GET_HOST_DATA_RDONLY(SELF%AIRD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSTAR)
#else
         !$acc exit data detach(SELF%WSTAR)
#endif
         CALL SELF%F_WSTAR%GET_HOST_DATA_RDONLY(SELF%WSTAR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CICOVER)
#else
         !$acc exit data detach(SELF%CICOVER)
#endif
         CALL SELF%F_CICOVER%GET_HOST_DATA_RDONLY(SELF%CICOVER)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CITHICK)
#else
         !$acc exit data detach(SELF%CITHICK)
#endif
         CALL SELF%F_CITHICK%GET_HOST_DATA_RDONLY(SELF%CITHICK)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%LKFR)
#else
         !$acc exit data detach(SELF%LKFR)
#endif
         CALL SELF%F_LKFR%GET_HOST_DATA_RDONLY(SELF%LKFR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%USTRA)
#else
         !$acc exit data detach(SELF%USTRA)
#endif
         CALL SELF%F_USTRA%GET_HOST_DATA_RDONLY(SELF%USTRA)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VSTRA)
#else
         !$acc exit data detach(SELF%VSTRA)
#endif
         CALL SELF%F_VSTRA%GET_HOST_DATA_RDONLY(SELF%VSTRA)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%UCUR)
#else
         !$acc exit data detach(SELF%UCUR)
#endif
         CALL SELF%F_UCUR%GET_HOST_DATA_RDONLY(SELF%UCUR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VCUR)
#else
         !$acc exit data detach(SELF%VCUR)
#endif
         CALL SELF%F_VCUR%GET_HOST_DATA_RDONLY(SELF%VCUR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSWAVE)
#else
         !$acc exit data detach(SELF%WSWAVE)
#endif
         CALL SELF%F_WSWAVE%GET_HOST_DATA_RDONLY(SELF%WSWAVE)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WDWAVE)
#else
         !$acc exit data detach(SELF%WDWAVE)
#endif
         CALL SELF%F_WDWAVE%GET_HOST_DATA_RDONLY(SELF%WDWAVE)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%UFRIC)
#else
         !$acc exit data detach(SELF%UFRIC)
#endif
         CALL SELF%F_UFRIC%GET_HOST_DATA_RDONLY(SELF%UFRIC)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUW)
#else
         !$acc exit data detach(SELF%TAUW)
#endif
         CALL SELF%F_TAUW%GET_HOST_DATA_RDONLY(SELF%TAUW)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUWDIR)
#else
         !$acc exit data detach(SELF%TAUWDIR)
#endif
         CALL SELF%F_TAUWDIR%GET_HOST_DATA_RDONLY(SELF%TAUWDIR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%Z0M)
#else
         !$acc exit data detach(SELF%Z0M)
#endif
         CALL SELF%F_Z0M%GET_HOST_DATA_RDONLY(SELF%Z0M)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%Z0B)
#else
         !$acc exit data detach(SELF%Z0B)
#endif
         CALL SELF%F_Z0B%GET_HOST_DATA_RDONLY(SELF%Z0B)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CHRNCK)
#else
         !$acc exit data detach(SELF%CHRNCK)
#endif
         CALL SELF%F_CHRNCK%GET_HOST_DATA_RDONLY(SELF%CHRNCK)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%XLON)
#else
         !$acc exit data detach(SELF%XLON)
#endif
         CALL SELF%F_XLON%GET_HOST_DATA_RDONLY(SELF%XLON)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%YLAT)
#else
         !$acc exit data detach(SELF%YLAT)
#endif
         CALL SELF%F_YLAT%GET_HOST_DATA_RDONLY(SELF%YLAT)
      ENDIF

   END SUBROUTINE FORCING_FIELDS_GET_HOST_DATA_RDONLY

   SUBROUTINE FORCING_FIELDS_GET_HOST_DATA_RDWR(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR, VCUR,&
       & WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT

      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%UWND)
#else
             !$acc exit data detach(SELF%UWND)
#endif
             CALL SELF%F_UWND%GET_HOST_DATA_RDWR(SELF%UWND)
         ENDIF
         IF(L_VWND)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VWND)
#else
             !$acc exit data detach(SELF%VWND)
#endif
             CALL SELF%F_VWND%GET_HOST_DATA_RDWR(SELF%VWND)
         ENDIF
         IF(L_AIRD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%AIRD)
#else
             !$acc exit data detach(SELF%AIRD)
#endif
             CALL SELF%F_AIRD%GET_HOST_DATA_RDWR(SELF%AIRD)
         ENDIF
         IF(L_WSTAR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSTAR)
#else
             !$acc exit data detach(SELF%WSTAR)
#endif
             CALL SELF%F_WSTAR%GET_HOST_DATA_RDWR(SELF%WSTAR)
         ENDIF
         IF(L_CICOVER)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CICOVER)
#else
             !$acc exit data detach(SELF%CICOVER)
#endif
             CALL SELF%F_CICOVER%GET_HOST_DATA_RDWR(SELF%CICOVER)
         ENDIF
         IF(L_CITHICK)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CITHICK)
#else
             !$acc exit data detach(SELF%CITHICK)
#endif
             CALL SELF%F_CITHICK%GET_HOST_DATA_RDWR(SELF%CITHICK)
         ENDIF
         IF(L_LKFR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%LKFR)
#else
             !$acc exit data detach(SELF%LKFR)
#endif
             CALL SELF%F_LKFR%GET_HOST_DATA_RDWR(SELF%LKFR)
         ENDIF
         IF(L_USTRA)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%USTRA)
#else
             !$acc exit data detach(SELF%USTRA)
#endif
             CALL SELF%F_USTRA%GET_HOST_DATA_RDWR(SELF%USTRA)
         ENDIF
         IF(L_VSTRA)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VSTRA)
#else
             !$acc exit data detach(SELF%VSTRA)
#endif
             CALL SELF%F_VSTRA%GET_HOST_DATA_RDWR(SELF%VSTRA)
         ENDIF
         IF(L_UCUR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%UCUR)
#else
             !$acc exit data detach(SELF%UCUR)
#endif
             CALL SELF%F_UCUR%GET_HOST_DATA_RDWR(SELF%UCUR)
         ENDIF
         IF(L_VCUR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VCUR)
#else
             !$acc exit data detach(SELF%VCUR)
#endif
             CALL SELF%F_VCUR%GET_HOST_DATA_RDWR(SELF%VCUR)
         ENDIF
         IF(L_WSWAVE)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSWAVE)
#else
             !$acc exit data detach(SELF%WSWAVE)
#endif
             CALL SELF%F_WSWAVE%GET_HOST_DATA_RDWR(SELF%WSWAVE)
         ENDIF
         IF(L_WDWAVE)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WDWAVE)
#else
             !$acc exit data detach(SELF%WDWAVE)
#endif
             CALL SELF%F_WDWAVE%GET_HOST_DATA_RDWR(SELF%WDWAVE)
         ENDIF
         IF(L_UFRIC)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%UFRIC)
#else
             !$acc exit data detach(SELF%UFRIC)
#endif
             CALL SELF%F_UFRIC%GET_HOST_DATA_RDWR(SELF%UFRIC)
         ENDIF
         IF(L_TAUW)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUW)
#else
             !$acc exit data detach(SELF%TAUW)
#endif
             CALL SELF%F_TAUW%GET_HOST_DATA_RDWR(SELF%TAUW)
         ENDIF
         IF(L_TAUWDIR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUWDIR)
#else
             !$acc exit data detach(SELF%TAUWDIR)
#endif
             CALL SELF%F_TAUWDIR%GET_HOST_DATA_RDWR(SELF%TAUWDIR)
         ENDIF
         IF(L_Z0M)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%Z0M)
#else
             !$acc exit data detach(SELF%Z0M)
#endif
             CALL SELF%F_Z0M%GET_HOST_DATA_RDWR(SELF%Z0M)
         ENDIF
         IF(L_Z0B)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%Z0B)
#else
             !$acc exit data detach(SELF%Z0B)
#endif
             CALL SELF%F_Z0B%GET_HOST_DATA_RDWR(SELF%Z0B)
         ENDIF
         IF(L_CHRNCK)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CHRNCK)
#else
             !$acc exit data detach(SELF%CHRNCK)
#endif
             CALL SELF%F_CHRNCK%GET_HOST_DATA_RDWR(SELF%CHRNCK)
         ENDIF
         IF(L_XLON)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%XLON)
#else
             !$acc exit data detach(SELF%XLON)
#endif
             CALL SELF%F_XLON%GET_HOST_DATA_RDWR(SELF%XLON)
         ENDIF
         IF(L_YLAT)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%YLAT)
#else
             !$acc exit data detach(SELF%YLAT)
#endif
             CALL SELF%F_YLAT%GET_HOST_DATA_RDWR(SELF%YLAT)
         ENDIF
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%UWND)
#else
         !$acc exit data detach(SELF%UWND)
#endif
         CALL SELF%F_UWND%GET_HOST_DATA_RDWR(SELF%UWND)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VWND)
#else
         !$acc exit data detach(SELF%VWND)
#endif
         CALL SELF%F_VWND%GET_HOST_DATA_RDWR(SELF%VWND)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%AIRD)
#else
         !$acc exit data detach(SELF%AIRD)
#endif
         CALL SELF%F_AIRD%GET_HOST_DATA_RDWR(SELF%AIRD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSTAR)
#else
         !$acc exit data detach(SELF%WSTAR)
#endif
         CALL SELF%F_WSTAR%GET_HOST_DATA_RDWR(SELF%WSTAR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CICOVER)
#else
         !$acc exit data detach(SELF%CICOVER)
#endif
         CALL SELF%F_CICOVER%GET_HOST_DATA_RDWR(SELF%CICOVER)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CITHICK)
#else
         !$acc exit data detach(SELF%CITHICK)
#endif
         CALL SELF%F_CITHICK%GET_HOST_DATA_RDWR(SELF%CITHICK)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%LKFR)
#else
         !$acc exit data detach(SELF%LKFR)
#endif
         CALL SELF%F_LKFR%GET_HOST_DATA_RDWR(SELF%LKFR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%USTRA)
#else
         !$acc exit data detach(SELF%USTRA)
#endif
         CALL SELF%F_USTRA%GET_HOST_DATA_RDWR(SELF%USTRA)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VSTRA)
#else
         !$acc exit data detach(SELF%VSTRA)
#endif
         CALL SELF%F_VSTRA%GET_HOST_DATA_RDWR(SELF%VSTRA)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%UCUR)
#else
         !$acc exit data detach(SELF%UCUR)
#endif
         CALL SELF%F_UCUR%GET_HOST_DATA_RDWR(SELF%UCUR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VCUR)
#else
         !$acc exit data detach(SELF%VCUR)
#endif
         CALL SELF%F_VCUR%GET_HOST_DATA_RDWR(SELF%VCUR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSWAVE)
#else
         !$acc exit data detach(SELF%WSWAVE)
#endif
         CALL SELF%F_WSWAVE%GET_HOST_DATA_RDWR(SELF%WSWAVE)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WDWAVE)
#else
         !$acc exit data detach(SELF%WDWAVE)
#endif
         CALL SELF%F_WDWAVE%GET_HOST_DATA_RDWR(SELF%WDWAVE)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%UFRIC)
#else
         !$acc exit data detach(SELF%UFRIC)
#endif
         CALL SELF%F_UFRIC%GET_HOST_DATA_RDWR(SELF%UFRIC)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUW)
#else
         !$acc exit data detach(SELF%TAUW)
#endif
         CALL SELF%F_TAUW%GET_HOST_DATA_RDWR(SELF%TAUW)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUWDIR)
#else
         !$acc exit data detach(SELF%TAUWDIR)
#endif
         CALL SELF%F_TAUWDIR%GET_HOST_DATA_RDWR(SELF%TAUWDIR)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%Z0M)
#else
         !$acc exit data detach(SELF%Z0M)
#endif
         CALL SELF%F_Z0M%GET_HOST_DATA_RDWR(SELF%Z0M)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%Z0B)
#else
         !$acc exit data detach(SELF%Z0B)
#endif
         CALL SELF%F_Z0B%GET_HOST_DATA_RDWR(SELF%Z0B)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CHRNCK)
#else
         !$acc exit data detach(SELF%CHRNCK)
#endif
         CALL SELF%F_CHRNCK%GET_HOST_DATA_RDWR(SELF%CHRNCK)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%XLON)
#else
         !$acc exit data detach(SELF%XLON)
#endif
         CALL SELF%F_XLON%GET_HOST_DATA_RDWR(SELF%XLON)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%YLAT)
#else
         !$acc exit data detach(SELF%YLAT)
#endif
         CALL SELF%F_YLAT%GET_HOST_DATA_RDWR(SELF%YLAT)
      ENDIF

   END SUBROUTINE FORCING_FIELDS_GET_HOST_DATA_RDWR

   SUBROUTINE FORCING_FIELDS_SYNC_HOST_RDONLY(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR, VCUR,&
       & WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT, QUEUE)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
             CALL SELF%F_UWND%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VWND)THEN
             CALL SELF%F_VWND%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_AIRD)THEN
             CALL SELF%F_AIRD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WSTAR)THEN
             CALL SELF%F_WSTAR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CICOVER)THEN
             CALL SELF%F_CICOVER%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CITHICK)THEN
             CALL SELF%F_CITHICK%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_LKFR)THEN
             CALL SELF%F_LKFR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_USTRA)THEN
             CALL SELF%F_USTRA%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTRA)THEN
             CALL SELF%F_VSTRA%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WSWAVE)THEN
             CALL SELF%F_WSWAVE%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WDWAVE)THEN
             CALL SELF%F_WDWAVE%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_UFRIC)THEN
             CALL SELF%F_UFRIC%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUW)THEN
             CALL SELF%F_TAUW%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUWDIR)THEN
             CALL SELF%F_TAUWDIR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0M)THEN
             CALL SELF%F_Z0M%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0B)THEN
             CALL SELF%F_Z0B%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CHRNCK)THEN
             CALL SELF%F_CHRNCK%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_XLON)THEN
             CALL SELF%F_XLON%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_YLAT)THEN
             CALL SELF%F_YLAT%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_UWND%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VWND%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_AIRD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WSTAR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CICOVER%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CITHICK%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_LKFR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_USTRA%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VSTRA%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WSWAVE%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WDWAVE%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_UFRIC%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUW%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUWDIR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_Z0M%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_Z0B%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CHRNCK%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_XLON%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_YLAT%SYNC_HOST_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE FORCING_FIELDS_SYNC_HOST_RDONLY

   SUBROUTINE FORCING_FIELDS_SYNC_HOST_RDWR(SELF, UWND, VWND, AIRD, WSTAR, CICOVER, CITHICK, LKFR, USTRA, VSTRA, UCUR, VCUR,&
       & WSWAVE, WDWAVE, UFRIC, TAUW, TAUWDIR, Z0M, Z0B, CHRNCK, XLON, YLAT, QUEUE)
      TYPE(FORCING_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: UWND
      LOGICAL :: L_UWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VWND
      LOGICAL :: L_VWND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: AIRD
      LOGICAL :: L_AIRD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSTAR
      LOGICAL :: L_WSTAR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CICOVER
      LOGICAL :: L_CICOVER = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CITHICK
      LOGICAL :: L_CITHICK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: LKFR
      LOGICAL :: L_LKFR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTRA
      LOGICAL :: L_USTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTRA
      LOGICAL :: L_VSTRA = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSWAVE
      LOGICAL :: L_WSWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WDWAVE
      LOGICAL :: L_WDWAVE = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UFRIC
      LOGICAL :: L_UFRIC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUW
      LOGICAL :: L_TAUW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUWDIR
      LOGICAL :: L_TAUWDIR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0M
      LOGICAL :: L_Z0M = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: Z0B
      LOGICAL :: L_Z0B = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CHRNCK
      LOGICAL :: L_CHRNCK = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLON
      LOGICAL :: L_XLON = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: YLAT
      LOGICAL :: L_YLAT = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(UWND)) L_UWND = UWND
     IF(PRESENT(VWND)) L_VWND = VWND
     IF(PRESENT(AIRD)) L_AIRD = AIRD
     IF(PRESENT(WSTAR)) L_WSTAR = WSTAR
     IF(PRESENT(CICOVER)) L_CICOVER = CICOVER
     IF(PRESENT(CITHICK)) L_CITHICK = CITHICK
     IF(PRESENT(LKFR)) L_LKFR = LKFR
     IF(PRESENT(USTRA)) L_USTRA = USTRA
     IF(PRESENT(VSTRA)) L_VSTRA = VSTRA
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR
     IF(PRESENT(WSWAVE)) L_WSWAVE = WSWAVE
     IF(PRESENT(WDWAVE)) L_WDWAVE = WDWAVE
     IF(PRESENT(UFRIC)) L_UFRIC = UFRIC
     IF(PRESENT(TAUW)) L_TAUW = TAUW
     IF(PRESENT(TAUWDIR)) L_TAUWDIR = TAUWDIR
     IF(PRESENT(Z0M)) L_Z0M = Z0M
     IF(PRESENT(Z0B)) L_Z0B = Z0B
     IF(PRESENT(CHRNCK)) L_CHRNCK = CHRNCK
     IF(PRESENT(XLON)) L_XLON = XLON
     IF(PRESENT(YLAT)) L_YLAT = YLAT


      IF(L_UWND.OR.L_VWND.OR.L_AIRD.OR.L_WSTAR.OR.L_CICOVER.OR.L_CITHICK.OR.L_LKFR.OR.L_USTRA.OR.L_VSTRA.OR.L_UCUR.OR.L_VCUR.OR.L_W&
          &SWAVE.OR.L_WDWAVE.OR.L_UFRIC.OR.L_TAUW.OR.L_TAUWDIR.OR.L_Z0M.OR.L_Z0B.OR.L_CHRNCK.OR.L_XLON.OR.L_YLAT)THEN
         !... copy only selected members
         IF(L_UWND)THEN
             CALL SELF%F_UWND%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VWND)THEN
             CALL SELF%F_VWND%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_AIRD)THEN
             CALL SELF%F_AIRD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WSTAR)THEN
             CALL SELF%F_WSTAR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CICOVER)THEN
             CALL SELF%F_CICOVER%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CITHICK)THEN
             CALL SELF%F_CITHICK%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_LKFR)THEN
             CALL SELF%F_LKFR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_USTRA)THEN
             CALL SELF%F_USTRA%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTRA)THEN
             CALL SELF%F_VSTRA%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WSWAVE)THEN
             CALL SELF%F_WSWAVE%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WDWAVE)THEN
             CALL SELF%F_WDWAVE%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_UFRIC)THEN
             CALL SELF%F_UFRIC%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUW)THEN
             CALL SELF%F_TAUW%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUWDIR)THEN
             CALL SELF%F_TAUWDIR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0M)THEN
             CALL SELF%F_Z0M%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_Z0B)THEN
             CALL SELF%F_Z0B%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CHRNCK)THEN
             CALL SELF%F_CHRNCK%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_XLON)THEN
             CALL SELF%F_XLON%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_YLAT)THEN
             CALL SELF%F_YLAT%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_UWND%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VWND%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_AIRD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WSTAR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CICOVER%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CITHICK%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_LKFR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_USTRA%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VSTRA%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WSWAVE%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WDWAVE%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_UFRIC%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUW%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUWDIR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_Z0M%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_Z0B%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CHRNCK%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_XLON%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_YLAT%SYNC_HOST_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE FORCING_FIELDS_SYNC_HOST_RDWR

   SUBROUTINE FORCING_FIELDS_DELETE_DEVICE_DATA(SELF)
      TYPE(FORCING_FIELDS) :: SELF

      CALL SELF%F_UWND%DELETE_DEVICE_DATA()
      CALL SELF%F_VWND%DELETE_DEVICE_DATA()
      CALL SELF%F_AIRD%DELETE_DEVICE_DATA()
      CALL SELF%F_WSTAR%DELETE_DEVICE_DATA()
      CALL SELF%F_CICOVER%DELETE_DEVICE_DATA()
      CALL SELF%F_CITHICK%DELETE_DEVICE_DATA()
      CALL SELF%F_LKFR%DELETE_DEVICE_DATA()
      CALL SELF%F_USTRA%DELETE_DEVICE_DATA()
      CALL SELF%F_VSTRA%DELETE_DEVICE_DATA()
      CALL SELF%F_UCUR%DELETE_DEVICE_DATA()
      CALL SELF%F_VCUR%DELETE_DEVICE_DATA()
      CALL SELF%F_WSWAVE%DELETE_DEVICE_DATA()
      CALL SELF%F_WDWAVE%DELETE_DEVICE_DATA()
      CALL SELF%F_UFRIC%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUW%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUWDIR%DELETE_DEVICE_DATA()
      CALL SELF%F_Z0M%DELETE_DEVICE_DATA()
      CALL SELF%F_Z0B%DELETE_DEVICE_DATA()
      CALL SELF%F_CHRNCK%DELETE_DEVICE_DATA()
      CALL SELF%F_XLON%DELETE_DEVICE_DATA()
      CALL SELF%F_YLAT%DELETE_DEVICE_DATA()

   END SUBROUTINE FORCING_FIELDS_DELETE_DEVICE_DATA

END MODULE FORCING_FIELDS_TYPE_MOD
