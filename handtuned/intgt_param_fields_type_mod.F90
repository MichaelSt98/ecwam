MODULE INTGT_PARAM_FIELDS_TYPE_MOD

   USE PARKIND_WAVE, ONLY : JWRB, JWIM, JWRO
   USE FIELD_MODULE, ONLY : FIELD_3RB, FIELD_2IM, FIELD_2RB, FIELD_3IM, FIELD_4RB, FIELD_4IM, FIELD_2RD, &
   &                        FIELD_2RM, FIELD_1IM, FIELD_1RB
   USE FIELD_FACTORY_MODULE, ONLY : FIELD_NEW, FIELD_DELETE
   IMPLICIT NONE

   PRIVATE

   INTERFACE ALLOC
      MODULE PROCEDURE INTGT_PARAM_FIELDS_ALLOC
   END INTERFACE ALLOC
   INTERFACE DEALLOC
      MODULE PROCEDURE INTGT_PARAM_FIELDS_DEALLOC
   END INTERFACE DEALLOC
   INTERFACE SYNC_DEVICE_RDWR
      MODULE PROCEDURE INTGT_PARAM_FIELDS_SYNC_DEVICE_RDWR
   END INTERFACE SYNC_DEVICE_RDWR
   INTERFACE SYNC_DEVICE_RDONLY
      MODULE PROCEDURE INTGT_PARAM_FIELDS_SYNC_DEVICE_RDONLY
   END INTERFACE SYNC_DEVICE_RDONLY
   INTERFACE SYNC_HOST_RDWR
      MODULE PROCEDURE INTGT_PARAM_FIELDS_SYNC_HOST_RDWR
   END INTERFACE SYNC_HOST_RDWR
   INTERFACE SYNC_HOST_RDONLY
      MODULE PROCEDURE INTGT_PARAM_FIELDS_SYNC_HOST_RDONLY
   END INTERFACE SYNC_HOST_RDONLY
   INTERFACE GET_DEVICE_DATA_RDWR
      MODULE PROCEDURE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_RDWR
   END INTERFACE GET_DEVICE_DATA_RDWR
   INTERFACE GET_DEVICE_DATA_RDONLY
      MODULE PROCEDURE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_RDONLY
   END INTERFACE GET_DEVICE_DATA_RDONLY
   INTERFACE GET_DEVICE_DATA_WRONLY
      MODULE PROCEDURE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_WRONLY
   END INTERFACE GET_DEVICE_DATA_WRONLY
   INTERFACE GET_HOST_DATA_RDWR
      MODULE PROCEDURE INTGT_PARAM_FIELDS_GET_HOST_DATA_RDWR
   END INTERFACE GET_HOST_DATA_RDWR
   INTERFACE GET_HOST_DATA_RDONLY
      MODULE PROCEDURE INTGT_PARAM_FIELDS_GET_HOST_DATA_RDONLY
   END INTERFACE GET_HOST_DATA_RDONLY
   INTERFACE DELETE_DEVICE_DATA
      MODULE PROCEDURE INTGT_PARAM_FIELDS_DELETE_DEVICE_DATA
   END INTERFACE DELETE_DEVICE_DATA

   TYPE INTGT_PARAM_FIELDS
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSEMEAN => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSEMEAN(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSEMEAN => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSEMEAN(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSFMEAN => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSFMEAN(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_WSFMEAN => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: WSFMEAN(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_USTOKES => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: USTOKES(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_USTOKES => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: USTOKES(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VSTOKES => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VSTOKES(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_VSTOKES => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: VSTOKES(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_PHIEPS => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIEPS(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_PHIEPS => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIEPS(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_PHIOCD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIOCD(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_PHIOCD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIOCD(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_PHIAW => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIAW(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_PHIAW => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: PHIAW(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUOC => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOC(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUOC => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOC(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUXD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUXD(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUXD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUXD(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUYD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUYD(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUYD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUYD(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUOCXD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCXD(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUOCXD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCXD(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUOCYD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCYD(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUOCYD => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUOCYD(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUICX => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUICX(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUICX => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUICX(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUICY => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUICY(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_TAUICY => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: TAUICY(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_STRNMS => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: STRNMS(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_STRNMS => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: STRNMS(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_ALTWH => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: ALTWH(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_ALTWH => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: ALTWH(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CALTWH => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CALTWH(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_CALTWH => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: CALTWH(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_RALTCOR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: RALTCOR(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_RALTCOR => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: RALTCOR(:,:) => NULL()
#endif
      LOGICAL :: LALLOC = .FALSE.
   END TYPE INTGT_PARAM_FIELDS

   PUBLIC :: INTGT_PARAM_FIELDS, ALLOC, DEALLOC
   PUBLIC :: SYNC_DEVICE_RDWR, SYNC_DEVICE_RDONLY, SYNC_HOST_RDWR, SYNC_HOST_RDONLY, GET_DEVICE_DATA_RDWR, GET_DEVICE_DATA_RDONLY,&
       & GET_DEVICE_DATA_WRONLY, GET_HOST_DATA_RDWR, GET_HOST_DATA_RDONLY, DELETE_DEVICE_DATA

   CONTAINS

   SUBROUTINE INTGT_PARAM_FIELDS_ALLOC(SELF, UBOUNDS, LBOUNDS)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      INTEGER(KIND=JWIM), INTENT(IN) :: UBOUNDS(2)
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: LBOUNDS(2)
      INTEGER(KIND=JWIM) :: LLBOUNDS(2)

      LLBOUNDS(:) = 1
      IF(PRESENT(LBOUNDS)) LLBOUNDS = LBOUNDS

      CALL FIELD_NEW(SELF%F_WSEMEAN, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_WSEMEAN%GET_HOST_DATA_RDWR(SELF%WSEMEAN)
      SELF%WSEMEAN(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_WSFMEAN, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_WSFMEAN%GET_HOST_DATA_RDWR(SELF%WSFMEAN)
      SELF%WSFMEAN(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_USTOKES, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_USTOKES%GET_HOST_DATA_RDWR(SELF%USTOKES)
      SELF%USTOKES(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_VSTOKES, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_VSTOKES%GET_HOST_DATA_RDWR(SELF%VSTOKES)
      SELF%VSTOKES(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_PHIEPS, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_PHIEPS%GET_HOST_DATA_RDWR(SELF%PHIEPS)
      SELF%PHIEPS(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_PHIOCD, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_PHIOCD%GET_HOST_DATA_RDWR(SELF%PHIOCD)
      SELF%PHIOCD(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_PHIAW, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_PHIAW%GET_HOST_DATA_RDWR(SELF%PHIAW)
      SELF%PHIAW(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUOC, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUOC%GET_HOST_DATA_RDWR(SELF%TAUOC)
      SELF%TAUOC(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUXD, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUXD%GET_HOST_DATA_RDWR(SELF%TAUXD)
      SELF%TAUXD(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUYD, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUYD%GET_HOST_DATA_RDWR(SELF%TAUYD)
      SELF%TAUYD(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUOCXD, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUOCXD%GET_HOST_DATA_RDWR(SELF%TAUOCXD)
      SELF%TAUOCXD(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUOCYD, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUOCYD%GET_HOST_DATA_RDWR(SELF%TAUOCYD)
      SELF%TAUOCYD(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUICX, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUICX%GET_HOST_DATA_RDWR(SELF%TAUICX)
      SELF%TAUICX(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_TAUICY, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_TAUICY%GET_HOST_DATA_RDWR(SELF%TAUICY)
      SELF%TAUICY(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_STRNMS, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_STRNMS%GET_HOST_DATA_RDWR(SELF%STRNMS)
      SELF%STRNMS(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_ALTWH, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_ALTWH%GET_HOST_DATA_RDWR(SELF%ALTWH)
      SELF%ALTWH(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_CALTWH, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_CALTWH%GET_HOST_DATA_RDWR(SELF%CALTWH)
      SELF%CALTWH(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_RALTCOR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_RALTCOR%GET_HOST_DATA_RDWR(SELF%RALTCOR)
      SELF%RALTCOR(:,:) = 0_JWRB

      SELF%LALLOC = .TRUE.
#ifdef OMPGPU
      !$omp target enter data map(to:SELF)
#else
      !$acc enter data copyin(SELF)
#endif

   END SUBROUTINE INTGT_PARAM_FIELDS_ALLOC

   SUBROUTINE INTGT_PARAM_FIELDS_DEALLOC(SELF)
      TYPE(INTGT_PARAM_FIELDS) :: SELF

      NULLIFY(SELF%WSEMEAN)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%WSEMEAN)
#else
      !$acc exit data detach(SELF%WSEMEAN)
#endif
      CALL FIELD_DELETE(SELF%F_WSEMEAN)
      NULLIFY(SELF%F_WSEMEAN)
      NULLIFY(SELF%WSFMEAN)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%WSFMEAN)
#else
      !$acc exit data detach(SELF%WSFMEAN)
#endif
      CALL FIELD_DELETE(SELF%F_WSFMEAN)
      NULLIFY(SELF%F_WSFMEAN)
      NULLIFY(SELF%USTOKES)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%USTOKES)
#else
      !$acc exit data detach(SELF%USTOKES)
#endif
      CALL FIELD_DELETE(SELF%F_USTOKES)
      NULLIFY(SELF%F_USTOKES)
      NULLIFY(SELF%VSTOKES)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%VSTOKES)
#else
      !$acc exit data detach(SELF%VSTOKES)
#endif
      CALL FIELD_DELETE(SELF%F_VSTOKES)
      NULLIFY(SELF%F_VSTOKES)
      NULLIFY(SELF%PHIEPS)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%PHIEPS)
#else
      !$acc exit data detach(SELF%PHIEPS)
#endif
      CALL FIELD_DELETE(SELF%F_PHIEPS)
      NULLIFY(SELF%F_PHIEPS)
      NULLIFY(SELF%PHIOCD)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%PHIOCD)
#else
      !$acc exit data detach(SELF%PHIOCD)
#endif
      CALL FIELD_DELETE(SELF%F_PHIOCD)
      NULLIFY(SELF%F_PHIOCD)
      NULLIFY(SELF%PHIAW)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%PHIAW)
#else
      !$acc exit data detach(SELF%PHIAW)
#endif
      CALL FIELD_DELETE(SELF%F_PHIAW)
      NULLIFY(SELF%F_PHIAW)
      NULLIFY(SELF%TAUOC)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUOC)
#else
      !$acc exit data detach(SELF%TAUOC)
#endif
      CALL FIELD_DELETE(SELF%F_TAUOC)
      NULLIFY(SELF%F_TAUOC)
      NULLIFY(SELF%TAUXD)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUXD)
#else
      !$acc exit data detach(SELF%TAUXD)
#endif
      CALL FIELD_DELETE(SELF%F_TAUXD)
      NULLIFY(SELF%F_TAUXD)
      NULLIFY(SELF%TAUYD)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUYD)
#else
      !$acc exit data detach(SELF%TAUYD)
#endif
      CALL FIELD_DELETE(SELF%F_TAUYD)
      NULLIFY(SELF%F_TAUYD)
      NULLIFY(SELF%TAUOCXD)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUOCXD)
#else
      !$acc exit data detach(SELF%TAUOCXD)
#endif
      CALL FIELD_DELETE(SELF%F_TAUOCXD)
      NULLIFY(SELF%F_TAUOCXD)
      NULLIFY(SELF%TAUOCYD)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUOCYD)
#else
      !$acc exit data detach(SELF%TAUOCYD)
#endif
      CALL FIELD_DELETE(SELF%F_TAUOCYD)
      NULLIFY(SELF%F_TAUOCYD)
      NULLIFY(SELF%TAUICX)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUICX)
#else
      !$acc exit data detach(SELF%TAUICX)
#endif
      CALL FIELD_DELETE(SELF%F_TAUICX)
      NULLIFY(SELF%F_TAUICX)
      NULLIFY(SELF%TAUICY)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%TAUICY)
#else
      !$acc exit data detach(SELF%TAUICY)
#endif
      CALL FIELD_DELETE(SELF%F_TAUICY)
      NULLIFY(SELF%F_TAUICY)
      NULLIFY(SELF%STRNMS)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%STRNMS)
#else
      !$acc exit data detach(SELF%STRNMS)
#endif
      CALL FIELD_DELETE(SELF%F_STRNMS)
      NULLIFY(SELF%F_STRNMS)
      NULLIFY(SELF%ALTWH)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%ALTWH)
#else
      !$acc exit data detach(SELF%ALTWH)
#endif
      CALL FIELD_DELETE(SELF%F_ALTWH)
      NULLIFY(SELF%F_ALTWH)
      NULLIFY(SELF%CALTWH)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%CALTWH)
#else
      !$acc exit data detach(SELF%CALTWH)
#endif
      CALL FIELD_DELETE(SELF%F_CALTWH)
      NULLIFY(SELF%F_CALTWH)
      NULLIFY(SELF%RALTCOR)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%RALTCOR)
#else
      !$acc exit data detach(SELF%RALTCOR)
#endif
      CALL FIELD_DELETE(SELF%F_RALTCOR)
      NULLIFY(SELF%F_RALTCOR)

      SELF%LALLOC = .FALSE.
#ifdef OMPGPU
      !$omp target exit data map(delete:SELF)
#else
      !$acc exit data delete(SELF)
#endif

   END SUBROUTINE INTGT_PARAM_FIELDS_DEALLOC

   SUBROUTINE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_RDWR(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC,&
       & TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
             CALL SELF%F_WSEMEAN%GET_DEVICE_DATA_RDWR(SELF%WSEMEAN)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSEMEAN)
#else
             !$acc enter data attach(SELF%WSEMEAN)
#endif
         ENDIF
         IF(L_WSFMEAN)THEN
             CALL SELF%F_WSFMEAN%GET_DEVICE_DATA_RDWR(SELF%WSFMEAN)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSFMEAN)
#else
             !$acc enter data attach(SELF%WSFMEAN)
#endif
         ENDIF
         IF(L_USTOKES)THEN
             CALL SELF%F_USTOKES%GET_DEVICE_DATA_RDWR(SELF%USTOKES)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%USTOKES)
#else
             !$acc enter data attach(SELF%USTOKES)
#endif
         ENDIF
         IF(L_VSTOKES)THEN
             CALL SELF%F_VSTOKES%GET_DEVICE_DATA_RDWR(SELF%VSTOKES)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VSTOKES)
#else
             !$acc enter data attach(SELF%VSTOKES)
#endif
         ENDIF
         IF(L_PHIEPS)THEN
             CALL SELF%F_PHIEPS%GET_DEVICE_DATA_RDWR(SELF%PHIEPS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIEPS)
#else
             !$acc enter data attach(SELF%PHIEPS)
#endif
         ENDIF
         IF(L_PHIOCD)THEN
             CALL SELF%F_PHIOCD%GET_DEVICE_DATA_RDWR(SELF%PHIOCD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIOCD)
#else
             !$acc enter data attach(SELF%PHIOCD)
#endif
         ENDIF
         IF(L_PHIAW)THEN
             CALL SELF%F_PHIAW%GET_DEVICE_DATA_RDWR(SELF%PHIAW)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIAW)
#else
             !$acc enter data attach(SELF%PHIAW)
#endif
         ENDIF
         IF(L_TAUOC)THEN
             CALL SELF%F_TAUOC%GET_DEVICE_DATA_RDWR(SELF%TAUOC)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOC)
#else
             !$acc enter data attach(SELF%TAUOC)
#endif
         ENDIF
         IF(L_TAUXD)THEN
             CALL SELF%F_TAUXD%GET_DEVICE_DATA_RDWR(SELF%TAUXD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUXD)
#else
             !$acc enter data attach(SELF%TAUXD)
#endif
         ENDIF
         IF(L_TAUYD)THEN
             CALL SELF%F_TAUYD%GET_DEVICE_DATA_RDWR(SELF%TAUYD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUYD)
#else
             !$acc enter data attach(SELF%TAUYD)
#endif
         ENDIF
         IF(L_TAUOCXD)THEN
             CALL SELF%F_TAUOCXD%GET_DEVICE_DATA_RDWR(SELF%TAUOCXD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOCXD)
#else
             !$acc enter data attach(SELF%TAUOCXD)
#endif
         ENDIF
         IF(L_TAUOCYD)THEN
             CALL SELF%F_TAUOCYD%GET_DEVICE_DATA_RDWR(SELF%TAUOCYD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOCYD)
#else
             !$acc enter data attach(SELF%TAUOCYD)
#endif
         ENDIF
         IF(L_TAUICX)THEN
             CALL SELF%F_TAUICX%GET_DEVICE_DATA_RDWR(SELF%TAUICX)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUICX)
#else
             !$acc enter data attach(SELF%TAUICX)
#endif
         ENDIF
         IF(L_TAUICY)THEN
             CALL SELF%F_TAUICY%GET_DEVICE_DATA_RDWR(SELF%TAUICY)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUICY)
#else
             !$acc enter data attach(SELF%TAUICY)
#endif
         ENDIF
         IF(L_STRNMS)THEN
             CALL SELF%F_STRNMS%GET_DEVICE_DATA_RDWR(SELF%STRNMS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%STRNMS)
#else
             !$acc enter data attach(SELF%STRNMS)
#endif
         ENDIF
         IF(L_ALTWH)THEN
             CALL SELF%F_ALTWH%GET_DEVICE_DATA_RDWR(SELF%ALTWH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%ALTWH)
#else
             !$acc enter data attach(SELF%ALTWH)
#endif
         ENDIF
         IF(L_CALTWH)THEN
             CALL SELF%F_CALTWH%GET_DEVICE_DATA_RDWR(SELF%CALTWH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CALTWH)
#else
             !$acc enter data attach(SELF%CALTWH)
#endif
         ENDIF
         IF(L_RALTCOR)THEN
             CALL SELF%F_RALTCOR%GET_DEVICE_DATA_RDWR(SELF%RALTCOR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%RALTCOR)
#else
             !$acc enter data attach(SELF%RALTCOR)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_WSEMEAN%GET_DEVICE_DATA_RDWR(SELF%WSEMEAN)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSEMEAN)
#else
         !$acc enter data attach(SELF%WSEMEAN)
#endif
         CALL SELF%F_WSFMEAN%GET_DEVICE_DATA_RDWR(SELF%WSFMEAN)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSFMEAN)
#else
         !$acc enter data attach(SELF%WSFMEAN)
#endif
         CALL SELF%F_USTOKES%GET_DEVICE_DATA_RDWR(SELF%USTOKES)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%USTOKES)
#else
         !$acc enter data attach(SELF%USTOKES)
#endif
         CALL SELF%F_VSTOKES%GET_DEVICE_DATA_RDWR(SELF%VSTOKES)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VSTOKES)
#else
         !$acc enter data attach(SELF%VSTOKES)
#endif
         CALL SELF%F_PHIEPS%GET_DEVICE_DATA_RDWR(SELF%PHIEPS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIEPS)
#else
         !$acc enter data attach(SELF%PHIEPS)
#endif
         CALL SELF%F_PHIOCD%GET_DEVICE_DATA_RDWR(SELF%PHIOCD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIOCD)
#else
         !$acc enter data attach(SELF%PHIOCD)
#endif
         CALL SELF%F_PHIAW%GET_DEVICE_DATA_RDWR(SELF%PHIAW)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIAW)
#else
         !$acc enter data attach(SELF%PHIAW)
#endif
         CALL SELF%F_TAUOC%GET_DEVICE_DATA_RDWR(SELF%TAUOC)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOC)
#else
         !$acc enter data attach(SELF%TAUOC)
#endif
         CALL SELF%F_TAUXD%GET_DEVICE_DATA_RDWR(SELF%TAUXD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUXD)
#else
         !$acc enter data attach(SELF%TAUXD)
#endif
         CALL SELF%F_TAUYD%GET_DEVICE_DATA_RDWR(SELF%TAUYD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUYD)
#else
         !$acc enter data attach(SELF%TAUYD)
#endif
         CALL SELF%F_TAUOCXD%GET_DEVICE_DATA_RDWR(SELF%TAUOCXD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOCXD)
#else
         !$acc enter data attach(SELF%TAUOCXD)
#endif
         CALL SELF%F_TAUOCYD%GET_DEVICE_DATA_RDWR(SELF%TAUOCYD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOCYD)
#else
         !$acc enter data attach(SELF%TAUOCYD)
#endif
         CALL SELF%F_TAUICX%GET_DEVICE_DATA_RDWR(SELF%TAUICX)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUICX)
#else
         !$acc enter data attach(SELF%TAUICX)
#endif
         CALL SELF%F_TAUICY%GET_DEVICE_DATA_RDWR(SELF%TAUICY)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUICY)
#else
         !$acc enter data attach(SELF%TAUICY)
#endif
         CALL SELF%F_STRNMS%GET_DEVICE_DATA_RDWR(SELF%STRNMS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%STRNMS)
#else
         !$acc enter data attach(SELF%STRNMS)
#endif
         CALL SELF%F_ALTWH%GET_DEVICE_DATA_RDWR(SELF%ALTWH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%ALTWH)
#else
         !$acc enter data attach(SELF%ALTWH)
#endif
         CALL SELF%F_CALTWH%GET_DEVICE_DATA_RDWR(SELF%CALTWH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CALTWH)
#else
         !$acc enter data attach(SELF%CALTWH)
#endif
         CALL SELF%F_RALTCOR%GET_DEVICE_DATA_RDWR(SELF%RALTCOR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%RALTCOR)
#else
         !$acc enter data attach(SELF%RALTCOR)
#endif
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_RDWR

   SUBROUTINE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_WRONLY(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC,&
       & TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
             CALL SELF%F_WSEMEAN%GET_DEVICE_DATA_WRONLY(SELF%WSEMEAN)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSEMEAN)
#else
             !$acc enter data attach(SELF%WSEMEAN)
#endif
         ENDIF
         IF(L_WSFMEAN)THEN
             CALL SELF%F_WSFMEAN%GET_DEVICE_DATA_WRONLY(SELF%WSFMEAN)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSFMEAN)
#else
             !$acc enter data attach(SELF%WSFMEAN)
#endif
         ENDIF
         IF(L_USTOKES)THEN
             CALL SELF%F_USTOKES%GET_DEVICE_DATA_WRONLY(SELF%USTOKES)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%USTOKES)
#else
             !$acc enter data attach(SELF%USTOKES)
#endif
         ENDIF
         IF(L_VSTOKES)THEN
             CALL SELF%F_VSTOKES%GET_DEVICE_DATA_WRONLY(SELF%VSTOKES)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VSTOKES)
#else
             !$acc enter data attach(SELF%VSTOKES)
#endif
         ENDIF
         IF(L_PHIEPS)THEN
             CALL SELF%F_PHIEPS%GET_DEVICE_DATA_WRONLY(SELF%PHIEPS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIEPS)
#else
             !$acc enter data attach(SELF%PHIEPS)
#endif
         ENDIF
         IF(L_PHIOCD)THEN
             CALL SELF%F_PHIOCD%GET_DEVICE_DATA_WRONLY(SELF%PHIOCD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIOCD)
#else
             !$acc enter data attach(SELF%PHIOCD)
#endif
         ENDIF
         IF(L_PHIAW)THEN
             CALL SELF%F_PHIAW%GET_DEVICE_DATA_WRONLY(SELF%PHIAW)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIAW)
#else
             !$acc enter data attach(SELF%PHIAW)
#endif
         ENDIF
         IF(L_TAUOC)THEN
             CALL SELF%F_TAUOC%GET_DEVICE_DATA_WRONLY(SELF%TAUOC)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOC)
#else
             !$acc enter data attach(SELF%TAUOC)
#endif
         ENDIF
         IF(L_TAUXD)THEN
             CALL SELF%F_TAUXD%GET_DEVICE_DATA_WRONLY(SELF%TAUXD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUXD)
#else
             !$acc enter data attach(SELF%TAUXD)
#endif
         ENDIF
         IF(L_TAUYD)THEN
             CALL SELF%F_TAUYD%GET_DEVICE_DATA_WRONLY(SELF%TAUYD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUYD)
#else
             !$acc enter data attach(SELF%TAUYD)
#endif
         ENDIF
         IF(L_TAUOCXD)THEN
             CALL SELF%F_TAUOCXD%GET_DEVICE_DATA_WRONLY(SELF%TAUOCXD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOCXD)
#else
             !$acc enter data attach(SELF%TAUOCXD)
#endif
         ENDIF
         IF(L_TAUOCYD)THEN
             CALL SELF%F_TAUOCYD%GET_DEVICE_DATA_WRONLY(SELF%TAUOCYD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOCYD)
#else
             !$acc enter data attach(SELF%TAUOCYD)
#endif
         ENDIF
         IF(L_TAUICX)THEN
             CALL SELF%F_TAUICX%GET_DEVICE_DATA_WRONLY(SELF%TAUICX)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUICX)
#else
             !$acc enter data attach(SELF%TAUICX)
#endif
         ENDIF
         IF(L_TAUICY)THEN
             CALL SELF%F_TAUICY%GET_DEVICE_DATA_WRONLY(SELF%TAUICY)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUICY)
#else
             !$acc enter data attach(SELF%TAUICY)
#endif
         ENDIF
         IF(L_STRNMS)THEN
             CALL SELF%F_STRNMS%GET_DEVICE_DATA_WRONLY(SELF%STRNMS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%STRNMS)
#else
             !$acc enter data attach(SELF%STRNMS)
#endif
         ENDIF
         IF(L_ALTWH)THEN
             CALL SELF%F_ALTWH%GET_DEVICE_DATA_WRONLY(SELF%ALTWH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%ALTWH)
#else
             !$acc enter data attach(SELF%ALTWH)
#endif
         ENDIF
         IF(L_CALTWH)THEN
             CALL SELF%F_CALTWH%GET_DEVICE_DATA_WRONLY(SELF%CALTWH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CALTWH)
#else
             !$acc enter data attach(SELF%CALTWH)
#endif
         ENDIF
         IF(L_RALTCOR)THEN
             CALL SELF%F_RALTCOR%GET_DEVICE_DATA_WRONLY(SELF%RALTCOR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%RALTCOR)
#else
             !$acc enter data attach(SELF%RALTCOR)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_WSEMEAN%GET_DEVICE_DATA_WRONLY(SELF%WSEMEAN)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSEMEAN)
#else
         !$acc enter data attach(SELF%WSEMEAN)
#endif
         CALL SELF%F_WSFMEAN%GET_DEVICE_DATA_WRONLY(SELF%WSFMEAN)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSFMEAN)
#else
         !$acc enter data attach(SELF%WSFMEAN)
#endif
         CALL SELF%F_USTOKES%GET_DEVICE_DATA_WRONLY(SELF%USTOKES)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%USTOKES)
#else
         !$acc enter data attach(SELF%USTOKES)
#endif
         CALL SELF%F_VSTOKES%GET_DEVICE_DATA_WRONLY(SELF%VSTOKES)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VSTOKES)
#else
         !$acc enter data attach(SELF%VSTOKES)
#endif
         CALL SELF%F_PHIEPS%GET_DEVICE_DATA_WRONLY(SELF%PHIEPS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIEPS)
#else
         !$acc enter data attach(SELF%PHIEPS)
#endif
         CALL SELF%F_PHIOCD%GET_DEVICE_DATA_WRONLY(SELF%PHIOCD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIOCD)
#else
         !$acc enter data attach(SELF%PHIOCD)
#endif
         CALL SELF%F_PHIAW%GET_DEVICE_DATA_WRONLY(SELF%PHIAW)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIAW)
#else
         !$acc enter data attach(SELF%PHIAW)
#endif
         CALL SELF%F_TAUOC%GET_DEVICE_DATA_WRONLY(SELF%TAUOC)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOC)
#else
         !$acc enter data attach(SELF%TAUOC)
#endif
         CALL SELF%F_TAUXD%GET_DEVICE_DATA_WRONLY(SELF%TAUXD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUXD)
#else
         !$acc enter data attach(SELF%TAUXD)
#endif
         CALL SELF%F_TAUYD%GET_DEVICE_DATA_WRONLY(SELF%TAUYD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUYD)
#else
         !$acc enter data attach(SELF%TAUYD)
#endif
         CALL SELF%F_TAUOCXD%GET_DEVICE_DATA_WRONLY(SELF%TAUOCXD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOCXD)
#else
         !$acc enter data attach(SELF%TAUOCXD)
#endif
         CALL SELF%F_TAUOCYD%GET_DEVICE_DATA_WRONLY(SELF%TAUOCYD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOCYD)
#else
         !$acc enter data attach(SELF%TAUOCYD)
#endif
         CALL SELF%F_TAUICX%GET_DEVICE_DATA_WRONLY(SELF%TAUICX)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUICX)
#else
         !$acc enter data attach(SELF%TAUICX)
#endif
         CALL SELF%F_TAUICY%GET_DEVICE_DATA_WRONLY(SELF%TAUICY)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUICY)
#else
         !$acc enter data attach(SELF%TAUICY)
#endif
         CALL SELF%F_STRNMS%GET_DEVICE_DATA_WRONLY(SELF%STRNMS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%STRNMS)
#else
         !$acc enter data attach(SELF%STRNMS)
#endif
         CALL SELF%F_ALTWH%GET_DEVICE_DATA_WRONLY(SELF%ALTWH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%ALTWH)
#else
         !$acc enter data attach(SELF%ALTWH)
#endif
         CALL SELF%F_CALTWH%GET_DEVICE_DATA_WRONLY(SELF%CALTWH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CALTWH)
#else
         !$acc enter data attach(SELF%CALTWH)
#endif
         CALL SELF%F_RALTCOR%GET_DEVICE_DATA_WRONLY(SELF%RALTCOR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%RALTCOR)
#else
         !$acc enter data attach(SELF%RALTCOR)
#endif
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_WRONLY

   SUBROUTINE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_RDONLY(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC,&
       & TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
             CALL SELF%F_WSEMEAN%GET_DEVICE_DATA_RDONLY(SELF%WSEMEAN)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSEMEAN)
#else
             !$acc enter data attach(SELF%WSEMEAN)
#endif
         ENDIF
         IF(L_WSFMEAN)THEN
             CALL SELF%F_WSFMEAN%GET_DEVICE_DATA_RDONLY(SELF%WSFMEAN)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%WSFMEAN)
#else
             !$acc enter data attach(SELF%WSFMEAN)
#endif
         ENDIF
         IF(L_USTOKES)THEN
             CALL SELF%F_USTOKES%GET_DEVICE_DATA_RDONLY(SELF%USTOKES)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%USTOKES)
#else
             !$acc enter data attach(SELF%USTOKES)
#endif
         ENDIF
         IF(L_VSTOKES)THEN
             CALL SELF%F_VSTOKES%GET_DEVICE_DATA_RDONLY(SELF%VSTOKES)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%VSTOKES)
#else
             !$acc enter data attach(SELF%VSTOKES)
#endif
         ENDIF
         IF(L_PHIEPS)THEN
             CALL SELF%F_PHIEPS%GET_DEVICE_DATA_RDONLY(SELF%PHIEPS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIEPS)
#else
             !$acc enter data attach(SELF%PHIEPS)
#endif
         ENDIF
         IF(L_PHIOCD)THEN
             CALL SELF%F_PHIOCD%GET_DEVICE_DATA_RDONLY(SELF%PHIOCD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIOCD)
#else
             !$acc enter data attach(SELF%PHIOCD)
#endif
         ENDIF
         IF(L_PHIAW)THEN
             CALL SELF%F_PHIAW%GET_DEVICE_DATA_RDONLY(SELF%PHIAW)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%PHIAW)
#else
             !$acc enter data attach(SELF%PHIAW)
#endif
         ENDIF
         IF(L_TAUOC)THEN
             CALL SELF%F_TAUOC%GET_DEVICE_DATA_RDONLY(SELF%TAUOC)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOC)
#else
             !$acc enter data attach(SELF%TAUOC)
#endif
         ENDIF
         IF(L_TAUXD)THEN
             CALL SELF%F_TAUXD%GET_DEVICE_DATA_RDONLY(SELF%TAUXD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUXD)
#else
             !$acc enter data attach(SELF%TAUXD)
#endif
         ENDIF
         IF(L_TAUYD)THEN
             CALL SELF%F_TAUYD%GET_DEVICE_DATA_RDONLY(SELF%TAUYD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUYD)
#else
             !$acc enter data attach(SELF%TAUYD)
#endif
         ENDIF
         IF(L_TAUOCXD)THEN
             CALL SELF%F_TAUOCXD%GET_DEVICE_DATA_RDONLY(SELF%TAUOCXD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOCXD)
#else
             !$acc enter data attach(SELF%TAUOCXD)
#endif
         ENDIF
         IF(L_TAUOCYD)THEN
             CALL SELF%F_TAUOCYD%GET_DEVICE_DATA_RDONLY(SELF%TAUOCYD)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUOCYD)
#else
             !$acc enter data attach(SELF%TAUOCYD)
#endif
         ENDIF
         IF(L_TAUICX)THEN
             CALL SELF%F_TAUICX%GET_DEVICE_DATA_RDONLY(SELF%TAUICX)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUICX)
#else
             !$acc enter data attach(SELF%TAUICX)
#endif
         ENDIF
         IF(L_TAUICY)THEN
             CALL SELF%F_TAUICY%GET_DEVICE_DATA_RDONLY(SELF%TAUICY)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%TAUICY)
#else
             !$acc enter data attach(SELF%TAUICY)
#endif
         ENDIF
         IF(L_STRNMS)THEN
             CALL SELF%F_STRNMS%GET_DEVICE_DATA_RDONLY(SELF%STRNMS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%STRNMS)
#else
             !$acc enter data attach(SELF%STRNMS)
#endif
         ENDIF
         IF(L_ALTWH)THEN
             CALL SELF%F_ALTWH%GET_DEVICE_DATA_RDONLY(SELF%ALTWH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%ALTWH)
#else
             !$acc enter data attach(SELF%ALTWH)
#endif
         ENDIF
         IF(L_CALTWH)THEN
             CALL SELF%F_CALTWH%GET_DEVICE_DATA_RDONLY(SELF%CALTWH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%CALTWH)
#else
             !$acc enter data attach(SELF%CALTWH)
#endif
         ENDIF
         IF(L_RALTCOR)THEN
             CALL SELF%F_RALTCOR%GET_DEVICE_DATA_RDONLY(SELF%RALTCOR)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%RALTCOR)
#else
             !$acc enter data attach(SELF%RALTCOR)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_WSEMEAN%GET_DEVICE_DATA_RDONLY(SELF%WSEMEAN)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSEMEAN)
#else
         !$acc enter data attach(SELF%WSEMEAN)
#endif
         CALL SELF%F_WSFMEAN%GET_DEVICE_DATA_RDONLY(SELF%WSFMEAN)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%WSFMEAN)
#else
         !$acc enter data attach(SELF%WSFMEAN)
#endif
         CALL SELF%F_USTOKES%GET_DEVICE_DATA_RDONLY(SELF%USTOKES)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%USTOKES)
#else
         !$acc enter data attach(SELF%USTOKES)
#endif
         CALL SELF%F_VSTOKES%GET_DEVICE_DATA_RDONLY(SELF%VSTOKES)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%VSTOKES)
#else
         !$acc enter data attach(SELF%VSTOKES)
#endif
         CALL SELF%F_PHIEPS%GET_DEVICE_DATA_RDONLY(SELF%PHIEPS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIEPS)
#else
         !$acc enter data attach(SELF%PHIEPS)
#endif
         CALL SELF%F_PHIOCD%GET_DEVICE_DATA_RDONLY(SELF%PHIOCD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIOCD)
#else
         !$acc enter data attach(SELF%PHIOCD)
#endif
         CALL SELF%F_PHIAW%GET_DEVICE_DATA_RDONLY(SELF%PHIAW)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%PHIAW)
#else
         !$acc enter data attach(SELF%PHIAW)
#endif
         CALL SELF%F_TAUOC%GET_DEVICE_DATA_RDONLY(SELF%TAUOC)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOC)
#else
         !$acc enter data attach(SELF%TAUOC)
#endif
         CALL SELF%F_TAUXD%GET_DEVICE_DATA_RDONLY(SELF%TAUXD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUXD)
#else
         !$acc enter data attach(SELF%TAUXD)
#endif
         CALL SELF%F_TAUYD%GET_DEVICE_DATA_RDONLY(SELF%TAUYD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUYD)
#else
         !$acc enter data attach(SELF%TAUYD)
#endif
         CALL SELF%F_TAUOCXD%GET_DEVICE_DATA_RDONLY(SELF%TAUOCXD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOCXD)
#else
         !$acc enter data attach(SELF%TAUOCXD)
#endif
         CALL SELF%F_TAUOCYD%GET_DEVICE_DATA_RDONLY(SELF%TAUOCYD)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUOCYD)
#else
         !$acc enter data attach(SELF%TAUOCYD)
#endif
         CALL SELF%F_TAUICX%GET_DEVICE_DATA_RDONLY(SELF%TAUICX)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUICX)
#else
         !$acc enter data attach(SELF%TAUICX)
#endif
         CALL SELF%F_TAUICY%GET_DEVICE_DATA_RDONLY(SELF%TAUICY)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%TAUICY)
#else
         !$acc enter data attach(SELF%TAUICY)
#endif
         CALL SELF%F_STRNMS%GET_DEVICE_DATA_RDONLY(SELF%STRNMS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%STRNMS)
#else
         !$acc enter data attach(SELF%STRNMS)
#endif
         CALL SELF%F_ALTWH%GET_DEVICE_DATA_RDONLY(SELF%ALTWH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%ALTWH)
#else
         !$acc enter data attach(SELF%ALTWH)
#endif
         CALL SELF%F_CALTWH%GET_DEVICE_DATA_RDONLY(SELF%CALTWH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%CALTWH)
#else
         !$acc enter data attach(SELF%CALTWH)
#endif
         CALL SELF%F_RALTCOR%GET_DEVICE_DATA_RDONLY(SELF%RALTCOR)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%RALTCOR)
#else
         !$acc enter data attach(SELF%RALTCOR)
#endif
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_GET_DEVICE_DATA_RDONLY

   SUBROUTINE INTGT_PARAM_FIELDS_SYNC_DEVICE_RDWR(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC, TAUXD,&
       & TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR, QUEUE)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
             CALL SELF%F_WSEMEAN%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WSFMEAN)THEN
             CALL SELF%F_WSFMEAN%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_USTOKES)THEN
             CALL SELF%F_USTOKES%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTOKES)THEN
             CALL SELF%F_VSTOKES%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIEPS)THEN
             CALL SELF%F_PHIEPS%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIOCD)THEN
             CALL SELF%F_PHIOCD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIAW)THEN
             CALL SELF%F_PHIAW%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOC)THEN
             CALL SELF%F_TAUOC%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUXD)THEN
             CALL SELF%F_TAUXD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUYD)THEN
             CALL SELF%F_TAUYD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCXD)THEN
             CALL SELF%F_TAUOCXD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCYD)THEN
             CALL SELF%F_TAUOCYD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICX)THEN
             CALL SELF%F_TAUICX%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICY)THEN
             CALL SELF%F_TAUICY%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_STRNMS)THEN
             CALL SELF%F_STRNMS%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_ALTWH)THEN
             CALL SELF%F_ALTWH%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CALTWH)THEN
             CALL SELF%F_CALTWH%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_RALTCOR)THEN
             CALL SELF%F_RALTCOR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_WSEMEAN%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WSFMEAN%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_USTOKES%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VSTOKES%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_PHIEPS%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_PHIOCD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_PHIAW%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUOC%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUXD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUYD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUOCXD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUOCYD%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUICX%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUICY%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_STRNMS%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_ALTWH%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CALTWH%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_RALTCOR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_SYNC_DEVICE_RDWR

   SUBROUTINE INTGT_PARAM_FIELDS_SYNC_DEVICE_RDONLY(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC, TAUXD,&
       & TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR, QUEUE)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
             CALL SELF%F_WSEMEAN%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WSFMEAN)THEN
             CALL SELF%F_WSFMEAN%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_USTOKES)THEN
             CALL SELF%F_USTOKES%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTOKES)THEN
             CALL SELF%F_VSTOKES%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIEPS)THEN
             CALL SELF%F_PHIEPS%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIOCD)THEN
             CALL SELF%F_PHIOCD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIAW)THEN
             CALL SELF%F_PHIAW%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOC)THEN
             CALL SELF%F_TAUOC%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUXD)THEN
             CALL SELF%F_TAUXD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUYD)THEN
             CALL SELF%F_TAUYD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCXD)THEN
             CALL SELF%F_TAUOCXD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCYD)THEN
             CALL SELF%F_TAUOCYD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICX)THEN
             CALL SELF%F_TAUICX%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICY)THEN
             CALL SELF%F_TAUICY%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_STRNMS)THEN
             CALL SELF%F_STRNMS%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_ALTWH)THEN
             CALL SELF%F_ALTWH%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CALTWH)THEN
             CALL SELF%F_CALTWH%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_RALTCOR)THEN
             CALL SELF%F_RALTCOR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_WSEMEAN%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WSFMEAN%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_USTOKES%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VSTOKES%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_PHIEPS%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_PHIOCD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_PHIAW%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUOC%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUXD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUYD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUOCXD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUOCYD%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUICX%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUICY%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_STRNMS%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_ALTWH%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CALTWH%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_RALTCOR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_SYNC_DEVICE_RDONLY

   SUBROUTINE INTGT_PARAM_FIELDS_GET_HOST_DATA_RDONLY(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC,&
       & TAUXD, TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSEMEAN)
#else
             !$acc exit data detach(SELF%WSEMEAN)
#endif
             CALL SELF%F_WSEMEAN%GET_HOST_DATA_RDONLY(SELF%WSEMEAN)
         ENDIF
         IF(L_WSFMEAN)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSFMEAN)
#else
             !$acc exit data detach(SELF%WSFMEAN)
#endif
             CALL SELF%F_WSFMEAN%GET_HOST_DATA_RDONLY(SELF%WSFMEAN)
         ENDIF
         IF(L_USTOKES)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%USTOKES)
#else
             !$acc exit data detach(SELF%USTOKES)
#endif
             CALL SELF%F_USTOKES%GET_HOST_DATA_RDONLY(SELF%USTOKES)
         ENDIF
         IF(L_VSTOKES)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VSTOKES)
#else
             !$acc exit data detach(SELF%VSTOKES)
#endif
             CALL SELF%F_VSTOKES%GET_HOST_DATA_RDONLY(SELF%VSTOKES)
         ENDIF
         IF(L_PHIEPS)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%PHIEPS)
#else
             !$acc exit data detach(SELF%PHIEPS)
#endif
             CALL SELF%F_PHIEPS%GET_HOST_DATA_RDONLY(SELF%PHIEPS)
         ENDIF
         IF(L_PHIOCD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%PHIOCD)
#else
             !$acc exit data detach(SELF%PHIOCD)
#endif
             CALL SELF%F_PHIOCD%GET_HOST_DATA_RDONLY(SELF%PHIOCD)
         ENDIF
         IF(L_PHIAW)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%PHIAW)
#else
             !$acc exit data detach(SELF%PHIAW)
#endif
             CALL SELF%F_PHIAW%GET_HOST_DATA_RDONLY(SELF%PHIAW)
         ENDIF
         IF(L_TAUOC)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUOC)
#else
             !$acc exit data detach(SELF%TAUOC)
#endif
             CALL SELF%F_TAUOC%GET_HOST_DATA_RDONLY(SELF%TAUOC)
         ENDIF
         IF(L_TAUXD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUXD)
#else
             !$acc exit data detach(SELF%TAUXD)
#endif
             CALL SELF%F_TAUXD%GET_HOST_DATA_RDONLY(SELF%TAUXD)
         ENDIF
         IF(L_TAUYD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUYD)
#else
             !$acc exit data detach(SELF%TAUYD)
#endif
             CALL SELF%F_TAUYD%GET_HOST_DATA_RDONLY(SELF%TAUYD)
         ENDIF
         IF(L_TAUOCXD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUOCXD)
#else
             !$acc exit data detach(SELF%TAUOCXD)
#endif
             CALL SELF%F_TAUOCXD%GET_HOST_DATA_RDONLY(SELF%TAUOCXD)
         ENDIF
         IF(L_TAUOCYD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUOCYD)
#else
             !$acc exit data detach(SELF%TAUOCYD)
#endif
             CALL SELF%F_TAUOCYD%GET_HOST_DATA_RDONLY(SELF%TAUOCYD)
         ENDIF
         IF(L_TAUICX)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUICX)
#else
             !$acc exit data detach(SELF%TAUICX)
#endif
             CALL SELF%F_TAUICX%GET_HOST_DATA_RDONLY(SELF%TAUICX)
         ENDIF
         IF(L_TAUICY)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUICY)
#else
             !$acc exit data detach(SELF%TAUICY)
#endif
             CALL SELF%F_TAUICY%GET_HOST_DATA_RDONLY(SELF%TAUICY)
         ENDIF
         IF(L_STRNMS)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%STRNMS)
#else
             !$acc exit data detach(SELF%STRNMS)
#endif
             CALL SELF%F_STRNMS%GET_HOST_DATA_RDONLY(SELF%STRNMS)
         ENDIF
         IF(L_ALTWH)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%ALTWH)
#else
             !$acc exit data detach(SELF%ALTWH)
#endif
             CALL SELF%F_ALTWH%GET_HOST_DATA_RDONLY(SELF%ALTWH)
         ENDIF
         IF(L_CALTWH)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CALTWH)
#else
             !$acc exit data detach(SELF%CALTWH)
#endif
             CALL SELF%F_CALTWH%GET_HOST_DATA_RDONLY(SELF%CALTWH)
         ENDIF
         IF(L_RALTCOR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%RALTCOR)
#else
             !$acc exit data detach(SELF%RALTCOR)
#endif
             CALL SELF%F_RALTCOR%GET_HOST_DATA_RDONLY(SELF%RALTCOR)
         ENDIF
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSEMEAN)
#else
         !$acc exit data detach(SELF%WSEMEAN)
#endif
         CALL SELF%F_WSEMEAN%GET_HOST_DATA_RDONLY(SELF%WSEMEAN)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSFMEAN)
#else
         !$acc exit data detach(SELF%WSFMEAN)
#endif
         CALL SELF%F_WSFMEAN%GET_HOST_DATA_RDONLY(SELF%WSFMEAN)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%USTOKES)
#else
         !$acc exit data detach(SELF%USTOKES)
#endif
         CALL SELF%F_USTOKES%GET_HOST_DATA_RDONLY(SELF%USTOKES)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VSTOKES)
#else
         !$acc exit data detach(SELF%VSTOKES)
#endif
         CALL SELF%F_VSTOKES%GET_HOST_DATA_RDONLY(SELF%VSTOKES)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%PHIEPS)
#else
         !$acc exit data detach(SELF%PHIEPS)
#endif
         CALL SELF%F_PHIEPS%GET_HOST_DATA_RDONLY(SELF%PHIEPS)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%PHIOCD)
#else
         !$acc exit data detach(SELF%PHIOCD)
#endif
         CALL SELF%F_PHIOCD%GET_HOST_DATA_RDONLY(SELF%PHIOCD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%PHIAW)
#else
         !$acc exit data detach(SELF%PHIAW)
#endif
         CALL SELF%F_PHIAW%GET_HOST_DATA_RDONLY(SELF%PHIAW)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUOC)
#else
         !$acc exit data detach(SELF%TAUOC)
#endif
         CALL SELF%F_TAUOC%GET_HOST_DATA_RDONLY(SELF%TAUOC)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUXD)
#else
         !$acc exit data detach(SELF%TAUXD)
#endif
         CALL SELF%F_TAUXD%GET_HOST_DATA_RDONLY(SELF%TAUXD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUYD)
#else
         !$acc exit data detach(SELF%TAUYD)
#endif
         CALL SELF%F_TAUYD%GET_HOST_DATA_RDONLY(SELF%TAUYD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUOCXD)
#else
         !$acc exit data detach(SELF%TAUOCXD)
#endif
         CALL SELF%F_TAUOCXD%GET_HOST_DATA_RDONLY(SELF%TAUOCXD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUOCYD)
#else
         !$acc exit data detach(SELF%TAUOCYD)
#endif
         CALL SELF%F_TAUOCYD%GET_HOST_DATA_RDONLY(SELF%TAUOCYD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUICX)
#else
         !$acc exit data detach(SELF%TAUICX)
#endif
         CALL SELF%F_TAUICX%GET_HOST_DATA_RDONLY(SELF%TAUICX)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUICY)
#else
         !$acc exit data detach(SELF%TAUICY)
#endif
         CALL SELF%F_TAUICY%GET_HOST_DATA_RDONLY(SELF%TAUICY)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%STRNMS)
#else
         !$acc exit data detach(SELF%STRNMS)
#endif
         CALL SELF%F_STRNMS%GET_HOST_DATA_RDONLY(SELF%STRNMS)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%ALTWH)
#else
         !$acc exit data detach(SELF%ALTWH)
#endif
         CALL SELF%F_ALTWH%GET_HOST_DATA_RDONLY(SELF%ALTWH)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CALTWH)
#else
         !$acc exit data detach(SELF%CALTWH)
#endif
         CALL SELF%F_CALTWH%GET_HOST_DATA_RDONLY(SELF%CALTWH)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%RALTCOR)
#else
         !$acc exit data detach(SELF%RALTCOR)
#endif
         CALL SELF%F_RALTCOR%GET_HOST_DATA_RDONLY(SELF%RALTCOR)
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_GET_HOST_DATA_RDONLY

   SUBROUTINE INTGT_PARAM_FIELDS_GET_HOST_DATA_RDWR(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC, TAUXD,&
       & TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR

      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSEMEAN)
#else
             !$acc exit data detach(SELF%WSEMEAN)
#endif
             CALL SELF%F_WSEMEAN%GET_HOST_DATA_RDWR(SELF%WSEMEAN)
         ENDIF
         IF(L_WSFMEAN)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%WSFMEAN)
#else
             !$acc exit data detach(SELF%WSFMEAN)
#endif
             CALL SELF%F_WSFMEAN%GET_HOST_DATA_RDWR(SELF%WSFMEAN)
         ENDIF
         IF(L_USTOKES)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%USTOKES)
#else
             !$acc exit data detach(SELF%USTOKES)
#endif
             CALL SELF%F_USTOKES%GET_HOST_DATA_RDWR(SELF%USTOKES)
         ENDIF
         IF(L_VSTOKES)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%VSTOKES)
#else
             !$acc exit data detach(SELF%VSTOKES)
#endif
             CALL SELF%F_VSTOKES%GET_HOST_DATA_RDWR(SELF%VSTOKES)
         ENDIF
         IF(L_PHIEPS)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%PHIEPS)
#else
             !$acc exit data detach(SELF%PHIEPS)
#endif
             CALL SELF%F_PHIEPS%GET_HOST_DATA_RDWR(SELF%PHIEPS)
         ENDIF
         IF(L_PHIOCD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%PHIOCD)
#else
             !$acc exit data detach(SELF%PHIOCD)
#endif
             CALL SELF%F_PHIOCD%GET_HOST_DATA_RDWR(SELF%PHIOCD)
         ENDIF
         IF(L_PHIAW)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%PHIAW)
#else
             !$acc exit data detach(SELF%PHIAW)
#endif
             CALL SELF%F_PHIAW%GET_HOST_DATA_RDWR(SELF%PHIAW)
         ENDIF
         IF(L_TAUOC)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUOC)
#else
             !$acc exit data detach(SELF%TAUOC)
#endif
             CALL SELF%F_TAUOC%GET_HOST_DATA_RDWR(SELF%TAUOC)
         ENDIF
         IF(L_TAUXD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUXD)
#else
             !$acc exit data detach(SELF%TAUXD)
#endif
             CALL SELF%F_TAUXD%GET_HOST_DATA_RDWR(SELF%TAUXD)
         ENDIF
         IF(L_TAUYD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUYD)
#else
             !$acc exit data detach(SELF%TAUYD)
#endif
             CALL SELF%F_TAUYD%GET_HOST_DATA_RDWR(SELF%TAUYD)
         ENDIF
         IF(L_TAUOCXD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUOCXD)
#else
             !$acc exit data detach(SELF%TAUOCXD)
#endif
             CALL SELF%F_TAUOCXD%GET_HOST_DATA_RDWR(SELF%TAUOCXD)
         ENDIF
         IF(L_TAUOCYD)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUOCYD)
#else
             !$acc exit data detach(SELF%TAUOCYD)
#endif
             CALL SELF%F_TAUOCYD%GET_HOST_DATA_RDWR(SELF%TAUOCYD)
         ENDIF
         IF(L_TAUICX)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUICX)
#else
             !$acc exit data detach(SELF%TAUICX)
#endif
             CALL SELF%F_TAUICX%GET_HOST_DATA_RDWR(SELF%TAUICX)
         ENDIF
         IF(L_TAUICY)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%TAUICY)
#else
             !$acc exit data detach(SELF%TAUICY)
#endif
             CALL SELF%F_TAUICY%GET_HOST_DATA_RDWR(SELF%TAUICY)
         ENDIF
         IF(L_STRNMS)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%STRNMS)
#else
             !$acc exit data detach(SELF%STRNMS)
#endif
             CALL SELF%F_STRNMS%GET_HOST_DATA_RDWR(SELF%STRNMS)
         ENDIF
         IF(L_ALTWH)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%ALTWH)
#else
             !$acc exit data detach(SELF%ALTWH)
#endif
             CALL SELF%F_ALTWH%GET_HOST_DATA_RDWR(SELF%ALTWH)
         ENDIF
         IF(L_CALTWH)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%CALTWH)
#else
             !$acc exit data detach(SELF%CALTWH)
#endif
             CALL SELF%F_CALTWH%GET_HOST_DATA_RDWR(SELF%CALTWH)
         ENDIF
         IF(L_RALTCOR)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%RALTCOR)
#else
             !$acc exit data detach(SELF%RALTCOR)
#endif
             CALL SELF%F_RALTCOR%GET_HOST_DATA_RDWR(SELF%RALTCOR)
         ENDIF
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSEMEAN)
#else
         !$acc exit data detach(SELF%WSEMEAN)
#endif
         CALL SELF%F_WSEMEAN%GET_HOST_DATA_RDWR(SELF%WSEMEAN)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%WSFMEAN)
#else
         !$acc exit data detach(SELF%WSFMEAN)
#endif
         CALL SELF%F_WSFMEAN%GET_HOST_DATA_RDWR(SELF%WSFMEAN)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%USTOKES)
#else
         !$acc exit data detach(SELF%USTOKES)
#endif
         CALL SELF%F_USTOKES%GET_HOST_DATA_RDWR(SELF%USTOKES)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%VSTOKES)
#else
         !$acc exit data detach(SELF%VSTOKES)
#endif
         CALL SELF%F_VSTOKES%GET_HOST_DATA_RDWR(SELF%VSTOKES)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%PHIEPS)
#else
         !$acc exit data detach(SELF%PHIEPS)
#endif
         CALL SELF%F_PHIEPS%GET_HOST_DATA_RDWR(SELF%PHIEPS)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%PHIOCD)
#else
         !$acc exit data detach(SELF%PHIOCD)
#endif
         CALL SELF%F_PHIOCD%GET_HOST_DATA_RDWR(SELF%PHIOCD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%PHIAW)
#else
         !$acc exit data detach(SELF%PHIAW)
#endif
         CALL SELF%F_PHIAW%GET_HOST_DATA_RDWR(SELF%PHIAW)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUOC)
#else
         !$acc exit data detach(SELF%TAUOC)
#endif
         CALL SELF%F_TAUOC%GET_HOST_DATA_RDWR(SELF%TAUOC)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUXD)
#else
         !$acc exit data detach(SELF%TAUXD)
#endif
         CALL SELF%F_TAUXD%GET_HOST_DATA_RDWR(SELF%TAUXD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUYD)
#else
         !$acc exit data detach(SELF%TAUYD)
#endif
         CALL SELF%F_TAUYD%GET_HOST_DATA_RDWR(SELF%TAUYD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUOCXD)
#else
         !$acc exit data detach(SELF%TAUOCXD)
#endif
         CALL SELF%F_TAUOCXD%GET_HOST_DATA_RDWR(SELF%TAUOCXD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUOCYD)
#else
         !$acc exit data detach(SELF%TAUOCYD)
#endif
         CALL SELF%F_TAUOCYD%GET_HOST_DATA_RDWR(SELF%TAUOCYD)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUICX)
#else
         !$acc exit data detach(SELF%TAUICX)
#endif
         CALL SELF%F_TAUICX%GET_HOST_DATA_RDWR(SELF%TAUICX)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%TAUICY)
#else
         !$acc exit data detach(SELF%TAUICY)
#endif
         CALL SELF%F_TAUICY%GET_HOST_DATA_RDWR(SELF%TAUICY)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%STRNMS)
#else
         !$acc exit data detach(SELF%STRNMS)
#endif
         CALL SELF%F_STRNMS%GET_HOST_DATA_RDWR(SELF%STRNMS)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%ALTWH)
#else
         !$acc exit data detach(SELF%ALTWH)
#endif
         CALL SELF%F_ALTWH%GET_HOST_DATA_RDWR(SELF%ALTWH)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%CALTWH)
#else
         !$acc exit data detach(SELF%CALTWH)
#endif
         CALL SELF%F_CALTWH%GET_HOST_DATA_RDWR(SELF%CALTWH)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%RALTCOR)
#else
         !$acc exit data detach(SELF%RALTCOR)
#endif
         CALL SELF%F_RALTCOR%GET_HOST_DATA_RDWR(SELF%RALTCOR)
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_GET_HOST_DATA_RDWR

   SUBROUTINE INTGT_PARAM_FIELDS_SYNC_HOST_RDONLY(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC, TAUXD,&
       & TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR, QUEUE)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
             CALL SELF%F_WSEMEAN%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_WSFMEAN)THEN
             CALL SELF%F_WSFMEAN%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_USTOKES)THEN
             CALL SELF%F_USTOKES%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTOKES)THEN
             CALL SELF%F_VSTOKES%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIEPS)THEN
             CALL SELF%F_PHIEPS%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIOCD)THEN
             CALL SELF%F_PHIOCD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIAW)THEN
             CALL SELF%F_PHIAW%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOC)THEN
             CALL SELF%F_TAUOC%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUXD)THEN
             CALL SELF%F_TAUXD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUYD)THEN
             CALL SELF%F_TAUYD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCXD)THEN
             CALL SELF%F_TAUOCXD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCYD)THEN
             CALL SELF%F_TAUOCYD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICX)THEN
             CALL SELF%F_TAUICX%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICY)THEN
             CALL SELF%F_TAUICY%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_STRNMS)THEN
             CALL SELF%F_STRNMS%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_ALTWH)THEN
             CALL SELF%F_ALTWH%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_CALTWH)THEN
             CALL SELF%F_CALTWH%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_RALTCOR)THEN
             CALL SELF%F_RALTCOR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_WSEMEAN%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_WSFMEAN%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_USTOKES%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VSTOKES%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_PHIEPS%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_PHIOCD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_PHIAW%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUOC%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUXD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUYD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUOCXD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUOCYD%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUICX%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_TAUICY%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_STRNMS%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_ALTWH%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_CALTWH%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_RALTCOR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_SYNC_HOST_RDONLY

   SUBROUTINE INTGT_PARAM_FIELDS_SYNC_HOST_RDWR(SELF, WSEMEAN, WSFMEAN, USTOKES, VSTOKES, PHIEPS, PHIOCD, PHIAW, TAUOC, TAUXD,&
       & TAUYD, TAUOCXD, TAUOCYD, TAUICX, TAUICY, STRNMS, ALTWH, CALTWH, RALTCOR, QUEUE)
      TYPE(INTGT_PARAM_FIELDS) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: WSEMEAN
      LOGICAL :: L_WSEMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: WSFMEAN
      LOGICAL :: L_WSFMEAN = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: USTOKES
      LOGICAL :: L_USTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VSTOKES
      LOGICAL :: L_VSTOKES = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIEPS
      LOGICAL :: L_PHIEPS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIOCD
      LOGICAL :: L_PHIOCD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: PHIAW
      LOGICAL :: L_PHIAW = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOC
      LOGICAL :: L_TAUOC = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUXD
      LOGICAL :: L_TAUXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUYD
      LOGICAL :: L_TAUYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCXD
      LOGICAL :: L_TAUOCXD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUOCYD
      LOGICAL :: L_TAUOCYD = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICX
      LOGICAL :: L_TAUICX = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: TAUICY
      LOGICAL :: L_TAUICY = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: STRNMS
      LOGICAL :: L_STRNMS = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: ALTWH
      LOGICAL :: L_ALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: CALTWH
      LOGICAL :: L_CALTWH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: RALTCOR
      LOGICAL :: L_RALTCOR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(WSEMEAN)) L_WSEMEAN = WSEMEAN
     IF(PRESENT(WSFMEAN)) L_WSFMEAN = WSFMEAN
     IF(PRESENT(USTOKES)) L_USTOKES = USTOKES
     IF(PRESENT(VSTOKES)) L_VSTOKES = VSTOKES
     IF(PRESENT(PHIEPS)) L_PHIEPS = PHIEPS
     IF(PRESENT(PHIOCD)) L_PHIOCD = PHIOCD
     IF(PRESENT(PHIAW)) L_PHIAW = PHIAW
     IF(PRESENT(TAUOC)) L_TAUOC = TAUOC
     IF(PRESENT(TAUXD)) L_TAUXD = TAUXD
     IF(PRESENT(TAUYD)) L_TAUYD = TAUYD
     IF(PRESENT(TAUOCXD)) L_TAUOCXD = TAUOCXD
     IF(PRESENT(TAUOCYD)) L_TAUOCYD = TAUOCYD
     IF(PRESENT(TAUICX)) L_TAUICX = TAUICX
     IF(PRESENT(TAUICY)) L_TAUICY = TAUICY
     IF(PRESENT(STRNMS)) L_STRNMS = STRNMS
     IF(PRESENT(ALTWH)) L_ALTWH = ALTWH
     IF(PRESENT(CALTWH)) L_CALTWH = CALTWH
     IF(PRESENT(RALTCOR)) L_RALTCOR = RALTCOR


      IF(L_WSEMEAN.OR.L_WSFMEAN.OR.L_USTOKES.OR.L_VSTOKES.OR.L_PHIEPS.OR.L_PHIOCD.OR.L_PHIAW.OR.L_TAUOC.OR.L_TAUXD.OR.L_TAUYD.OR.L_&
          &TAUOCXD.OR.L_TAUOCYD.OR.L_TAUICX.OR.L_TAUICY.OR.L_STRNMS.OR.L_ALTWH.OR.L_CALTWH.OR.L_RALTCOR)THEN
         !... copy only selected members
         IF(L_WSEMEAN)THEN
             CALL SELF%F_WSEMEAN%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_WSFMEAN)THEN
             CALL SELF%F_WSFMEAN%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_USTOKES)THEN
             CALL SELF%F_USTOKES%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VSTOKES)THEN
             CALL SELF%F_VSTOKES%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIEPS)THEN
             CALL SELF%F_PHIEPS%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIOCD)THEN
             CALL SELF%F_PHIOCD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_PHIAW)THEN
             CALL SELF%F_PHIAW%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOC)THEN
             CALL SELF%F_TAUOC%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUXD)THEN
             CALL SELF%F_TAUXD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUYD)THEN
             CALL SELF%F_TAUYD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCXD)THEN
             CALL SELF%F_TAUOCXD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUOCYD)THEN
             CALL SELF%F_TAUOCYD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICX)THEN
             CALL SELF%F_TAUICX%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_TAUICY)THEN
             CALL SELF%F_TAUICY%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_STRNMS)THEN
             CALL SELF%F_STRNMS%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_ALTWH)THEN
             CALL SELF%F_ALTWH%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_CALTWH)THEN
             CALL SELF%F_CALTWH%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_RALTCOR)THEN
             CALL SELF%F_RALTCOR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_WSEMEAN%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_WSFMEAN%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_USTOKES%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VSTOKES%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_PHIEPS%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_PHIOCD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_PHIAW%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUOC%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUXD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUYD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUOCXD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUOCYD%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUICX%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_TAUICY%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_STRNMS%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_ALTWH%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_CALTWH%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_RALTCOR%SYNC_HOST_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE INTGT_PARAM_FIELDS_SYNC_HOST_RDWR

   SUBROUTINE INTGT_PARAM_FIELDS_DELETE_DEVICE_DATA(SELF)
      TYPE(INTGT_PARAM_FIELDS) :: SELF

      CALL SELF%F_WSEMEAN%DELETE_DEVICE_DATA()
      CALL SELF%F_WSFMEAN%DELETE_DEVICE_DATA()
      CALL SELF%F_USTOKES%DELETE_DEVICE_DATA()
      CALL SELF%F_VSTOKES%DELETE_DEVICE_DATA()
      CALL SELF%F_PHIEPS%DELETE_DEVICE_DATA()
      CALL SELF%F_PHIOCD%DELETE_DEVICE_DATA()
      CALL SELF%F_PHIAW%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUOC%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUXD%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUYD%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUOCXD%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUOCYD%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUICX%DELETE_DEVICE_DATA()
      CALL SELF%F_TAUICY%DELETE_DEVICE_DATA()
      CALL SELF%F_STRNMS%DELETE_DEVICE_DATA()
      CALL SELF%F_ALTWH%DELETE_DEVICE_DATA()
      CALL SELF%F_CALTWH%DELETE_DEVICE_DATA()
      CALL SELF%F_RALTCOR%DELETE_DEVICE_DATA()

   END SUBROUTINE INTGT_PARAM_FIELDS_DELETE_DEVICE_DATA

END MODULE INTGT_PARAM_FIELDS_TYPE_MOD
