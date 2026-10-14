MODULE ENVIRONMENT_TYPE_MOD

   USE PARKIND_WAVE, ONLY : JWRB, JWIM, JWRO
   USE FIELD_MODULE, ONLY : FIELD_3RB, FIELD_2IM, FIELD_2RB, FIELD_3IM, FIELD_4RB, FIELD_4IM, FIELD_2RD, &
   &                        FIELD_2RM, FIELD_1IM, FIELD_1RB
   USE FIELD_FACTORY_MODULE, ONLY : FIELD_NEW, FIELD_DELETE
   IMPLICIT NONE

   PRIVATE

   INTERFACE ALLOC
      MODULE PROCEDURE ENVIRONMENT_ALLOC
   END INTERFACE ALLOC
   INTERFACE DEALLOC
      MODULE PROCEDURE ENVIRONMENT_DEALLOC
   END INTERFACE DEALLOC
   INTERFACE SYNC_DEVICE_RDWR
      MODULE PROCEDURE ENVIRONMENT_SYNC_DEVICE_RDWR
   END INTERFACE SYNC_DEVICE_RDWR
   INTERFACE SYNC_DEVICE_RDONLY
      MODULE PROCEDURE ENVIRONMENT_SYNC_DEVICE_RDONLY
   END INTERFACE SYNC_DEVICE_RDONLY
   INTERFACE SYNC_HOST_RDWR
      MODULE PROCEDURE ENVIRONMENT_SYNC_HOST_RDWR
   END INTERFACE SYNC_HOST_RDWR
   INTERFACE SYNC_HOST_RDONLY
      MODULE PROCEDURE ENVIRONMENT_SYNC_HOST_RDONLY
   END INTERFACE SYNC_HOST_RDONLY
   INTERFACE GET_DEVICE_DATA_RDWR
      MODULE PROCEDURE ENVIRONMENT_GET_DEVICE_DATA_RDWR
   END INTERFACE GET_DEVICE_DATA_RDWR
   INTERFACE GET_DEVICE_DATA_RDONLY
      MODULE PROCEDURE ENVIRONMENT_GET_DEVICE_DATA_RDONLY
   END INTERFACE GET_DEVICE_DATA_RDONLY
   INTERFACE GET_DEVICE_DATA_WRONLY
      MODULE PROCEDURE ENVIRONMENT_GET_DEVICE_DATA_WRONLY
   END INTERFACE GET_DEVICE_DATA_WRONLY
   INTERFACE GET_HOST_DATA_RDWR
      MODULE PROCEDURE ENVIRONMENT_GET_HOST_DATA_RDWR
   END INTERFACE GET_HOST_DATA_RDWR
   INTERFACE GET_HOST_DATA_RDONLY
      MODULE PROCEDURE ENVIRONMENT_GET_HOST_DATA_RDONLY
   END INTERFACE GET_HOST_DATA_RDONLY
   INTERFACE DELETE_DEVICE_DATA
      MODULE PROCEDURE ENVIRONMENT_DELETE_DEVICE_DATA
   END INTERFACE DELETE_DEVICE_DATA

   TYPE ENVIRONMENT
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2IM), POINTER :: F_INDEP => NULL()
      INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: INDEP(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2IM), POINTER :: F_INDEP => NULL()
      INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: INDEP(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2IM), POINTER :: F_IODP => NULL()
      INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IODP(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2IM), POINTER :: F_IODP => NULL()
      INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IODP(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2IM), POINTER :: F_IOBND => NULL()
      INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IOBND(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2IM), POINTER :: F_IOBND => NULL()
      INTEGER(KIND=JWIM), POINTER, CONTIGUOUS :: IOBND(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_IBRMEM => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: IBRMEM(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_IBRMEM => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: IBRMEM(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_DELLAM1 => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: DELLAM1(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_DELLAM1 => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: DELLAM1(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_COSPHM1 => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: COSPHM1(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_COSPHM1 => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: COSPHM1(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_DEPTH => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: DEPTH(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_DEPTH => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: DEPTH(:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_EMAXDPT => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: EMAXDPT(:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_2RB), POINTER :: F_EMAXDPT => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: EMAXDPT(:,:) => NULL()
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
      LOGICAL :: LALLOC = .FALSE.
   END TYPE ENVIRONMENT

   PUBLIC :: ENVIRONMENT, ALLOC, DEALLOC
   PUBLIC :: SYNC_DEVICE_RDWR, SYNC_DEVICE_RDONLY, SYNC_HOST_RDWR, SYNC_HOST_RDONLY, GET_DEVICE_DATA_RDWR, GET_DEVICE_DATA_RDONLY,&
       & GET_DEVICE_DATA_WRONLY, GET_HOST_DATA_RDWR, GET_HOST_DATA_RDONLY, DELETE_DEVICE_DATA

   CONTAINS

   SUBROUTINE ENVIRONMENT_ALLOC(SELF, UBOUNDS, LBOUNDS)
      TYPE(ENVIRONMENT) :: SELF
      INTEGER(KIND=JWIM), INTENT(IN) :: UBOUNDS(2)
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: LBOUNDS(2)
      INTEGER(KIND=JWIM) :: LLBOUNDS(2)

      LLBOUNDS(:) = 1
      IF(PRESENT(LBOUNDS)) LLBOUNDS = LBOUNDS

      CALL FIELD_NEW(SELF%F_INDEP, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_INDEP%GET_HOST_DATA_RDWR(SELF%INDEP)
      SELF%INDEP(:,:) = 0_JWIM
      CALL FIELD_NEW(SELF%F_IODP, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_IODP%GET_HOST_DATA_RDWR(SELF%IODP)
      SELF%IODP(:,:) = 0_JWIM
      CALL FIELD_NEW(SELF%F_IOBND, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_IOBND%GET_HOST_DATA_RDWR(SELF%IOBND)
      SELF%IOBND(:,:) = 0_JWIM
      CALL FIELD_NEW(SELF%F_IBRMEM, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_IBRMEM%GET_HOST_DATA_RDWR(SELF%IBRMEM)
      SELF%IBRMEM(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_DELLAM1, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_DELLAM1%GET_HOST_DATA_RDWR(SELF%DELLAM1)
      SELF%DELLAM1(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_COSPHM1, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_COSPHM1%GET_HOST_DATA_RDWR(SELF%COSPHM1)
      SELF%COSPHM1(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_DEPTH, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_DEPTH%GET_HOST_DATA_RDWR(SELF%DEPTH)
      SELF%DEPTH(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_EMAXDPT, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_EMAXDPT%GET_HOST_DATA_RDWR(SELF%EMAXDPT)
      SELF%EMAXDPT(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_UCUR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_UCUR%GET_HOST_DATA_RDWR(SELF%UCUR)
      SELF%UCUR(:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_VCUR, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_VCUR%GET_HOST_DATA_RDWR(SELF%VCUR)
      SELF%VCUR(:,:) = 0_JWRB

      SELF%LALLOC = .TRUE.
#ifdef OMPGPU
      !$omp target enter data map(to:SELF)
#else
      !$acc enter data copyin(SELF)
#endif

   END SUBROUTINE ENVIRONMENT_ALLOC

   SUBROUTINE ENVIRONMENT_DEALLOC(SELF)
      TYPE(ENVIRONMENT) :: SELF

      NULLIFY(SELF%INDEP)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%INDEP)
#else
      !$acc exit data detach(SELF%INDEP)
#endif
      CALL FIELD_DELETE(SELF%F_INDEP)
      NULLIFY(SELF%F_INDEP)
      NULLIFY(SELF%IODP)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%IODP)
#else
      !$acc exit data detach(SELF%IODP)
#endif
      CALL FIELD_DELETE(SELF%F_IODP)
      NULLIFY(SELF%F_IODP)
      NULLIFY(SELF%IOBND)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%IOBND)
#else
      !$acc exit data detach(SELF%IOBND)
#endif
      CALL FIELD_DELETE(SELF%F_IOBND)
      NULLIFY(SELF%F_IOBND)
      NULLIFY(SELF%IBRMEM)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%IBRMEM)
#else
      !$acc exit data detach(SELF%IBRMEM)
#endif
      CALL FIELD_DELETE(SELF%F_IBRMEM)
      NULLIFY(SELF%F_IBRMEM)
      NULLIFY(SELF%DELLAM1)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%DELLAM1)
#else
      !$acc exit data detach(SELF%DELLAM1)
#endif
      CALL FIELD_DELETE(SELF%F_DELLAM1)
      NULLIFY(SELF%F_DELLAM1)
      NULLIFY(SELF%COSPHM1)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%COSPHM1)
#else
      !$acc exit data detach(SELF%COSPHM1)
#endif
      CALL FIELD_DELETE(SELF%F_COSPHM1)
      NULLIFY(SELF%F_COSPHM1)
      NULLIFY(SELF%DEPTH)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%DEPTH)
#else
      !$acc exit data detach(SELF%DEPTH)
#endif
      CALL FIELD_DELETE(SELF%F_DEPTH)
      NULLIFY(SELF%F_DEPTH)
      NULLIFY(SELF%EMAXDPT)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%EMAXDPT)
#else
      !$acc exit data detach(SELF%EMAXDPT)
#endif
      CALL FIELD_DELETE(SELF%F_EMAXDPT)
      NULLIFY(SELF%F_EMAXDPT)
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

      SELF%LALLOC = .FALSE.
#ifdef OMPGPU
      !$omp target exit data map(delete:SELF)
#else
      !$acc exit data delete(SELF)
#endif

   END SUBROUTINE ENVIRONMENT_DEALLOC

   SUBROUTINE ENVIRONMENT_GET_DEVICE_DATA_RDWR(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
             CALL SELF%F_INDEP%GET_DEVICE_DATA_RDWR(SELF%INDEP)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%INDEP)
#else
             !$acc enter data attach(SELF%INDEP)
#endif
         ENDIF
         IF(L_IODP)THEN
             CALL SELF%F_IODP%GET_DEVICE_DATA_RDWR(SELF%IODP)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IODP)
#else
             !$acc enter data attach(SELF%IODP)
#endif
         ENDIF
         IF(L_IOBND)THEN
             CALL SELF%F_IOBND%GET_DEVICE_DATA_RDWR(SELF%IOBND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IOBND)
#else
             !$acc enter data attach(SELF%IOBND)
#endif
         ENDIF
         IF(L_IBRMEM)THEN
             CALL SELF%F_IBRMEM%GET_DEVICE_DATA_RDWR(SELF%IBRMEM)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IBRMEM)
#else
             !$acc enter data attach(SELF%IBRMEM)
#endif
         ENDIF
         IF(L_DELLAM1)THEN
             CALL SELF%F_DELLAM1%GET_DEVICE_DATA_RDWR(SELF%DELLAM1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%DELLAM1)
#else
             !$acc enter data attach(SELF%DELLAM1)
#endif
         ENDIF
         IF(L_COSPHM1)THEN
             CALL SELF%F_COSPHM1%GET_DEVICE_DATA_RDWR(SELF%COSPHM1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%COSPHM1)
#else
             !$acc enter data attach(SELF%COSPHM1)
#endif
         ENDIF
         IF(L_DEPTH)THEN
             CALL SELF%F_DEPTH%GET_DEVICE_DATA_RDWR(SELF%DEPTH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%DEPTH)
#else
             !$acc enter data attach(SELF%DEPTH)
#endif
         ENDIF
         IF(L_EMAXDPT)THEN
             CALL SELF%F_EMAXDPT%GET_DEVICE_DATA_RDWR(SELF%EMAXDPT)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%EMAXDPT)
#else
             !$acc enter data attach(SELF%EMAXDPT)
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
      ELSE
         !... copy entire struct
         CALL SELF%F_INDEP%GET_DEVICE_DATA_RDWR(SELF%INDEP)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%INDEP)
#else
         !$acc enter data attach(SELF%INDEP)
#endif
         CALL SELF%F_IODP%GET_DEVICE_DATA_RDWR(SELF%IODP)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IODP)
#else
         !$acc enter data attach(SELF%IODP)
#endif
         CALL SELF%F_IOBND%GET_DEVICE_DATA_RDWR(SELF%IOBND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IOBND)
#else
         !$acc enter data attach(SELF%IOBND)
#endif
         CALL SELF%F_IBRMEM%GET_DEVICE_DATA_RDWR(SELF%IBRMEM)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IBRMEM)
#else
         !$acc enter data attach(SELF%IBRMEM)
#endif
         CALL SELF%F_DELLAM1%GET_DEVICE_DATA_RDWR(SELF%DELLAM1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%DELLAM1)
#else
         !$acc enter data attach(SELF%DELLAM1)
#endif
         CALL SELF%F_COSPHM1%GET_DEVICE_DATA_RDWR(SELF%COSPHM1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%COSPHM1)
#else
         !$acc enter data attach(SELF%COSPHM1)
#endif
         CALL SELF%F_DEPTH%GET_DEVICE_DATA_RDWR(SELF%DEPTH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%DEPTH)
#else
         !$acc enter data attach(SELF%DEPTH)
#endif
         CALL SELF%F_EMAXDPT%GET_DEVICE_DATA_RDWR(SELF%EMAXDPT)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%EMAXDPT)
#else
         !$acc enter data attach(SELF%EMAXDPT)
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
      ENDIF

   END SUBROUTINE ENVIRONMENT_GET_DEVICE_DATA_RDWR

   SUBROUTINE ENVIRONMENT_GET_DEVICE_DATA_WRONLY(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
             CALL SELF%F_INDEP%GET_DEVICE_DATA_WRONLY(SELF%INDEP)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%INDEP)
#else
             !$acc enter data attach(SELF%INDEP)
#endif
         ENDIF
         IF(L_IODP)THEN
             CALL SELF%F_IODP%GET_DEVICE_DATA_WRONLY(SELF%IODP)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IODP)
#else
             !$acc enter data attach(SELF%IODP)
#endif
         ENDIF
         IF(L_IOBND)THEN
             CALL SELF%F_IOBND%GET_DEVICE_DATA_WRONLY(SELF%IOBND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IOBND)
#else
             !$acc enter data attach(SELF%IOBND)
#endif
         ENDIF
         IF(L_IBRMEM)THEN
             CALL SELF%F_IBRMEM%GET_DEVICE_DATA_WRONLY(SELF%IBRMEM)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IBRMEM)
#else
             !$acc enter data attach(SELF%IBRMEM)
#endif
         ENDIF
         IF(L_DELLAM1)THEN
             CALL SELF%F_DELLAM1%GET_DEVICE_DATA_WRONLY(SELF%DELLAM1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%DELLAM1)
#else
             !$acc enter data attach(SELF%DELLAM1)
#endif
         ENDIF
         IF(L_COSPHM1)THEN
             CALL SELF%F_COSPHM1%GET_DEVICE_DATA_WRONLY(SELF%COSPHM1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%COSPHM1)
#else
             !$acc enter data attach(SELF%COSPHM1)
#endif
         ENDIF
         IF(L_DEPTH)THEN
             CALL SELF%F_DEPTH%GET_DEVICE_DATA_WRONLY(SELF%DEPTH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%DEPTH)
#else
             !$acc enter data attach(SELF%DEPTH)
#endif
         ENDIF
         IF(L_EMAXDPT)THEN
             CALL SELF%F_EMAXDPT%GET_DEVICE_DATA_WRONLY(SELF%EMAXDPT)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%EMAXDPT)
#else
             !$acc enter data attach(SELF%EMAXDPT)
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
      ELSE
         !... copy entire struct
         CALL SELF%F_INDEP%GET_DEVICE_DATA_WRONLY(SELF%INDEP)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%INDEP)
#else
         !$acc enter data attach(SELF%INDEP)
#endif
         CALL SELF%F_IODP%GET_DEVICE_DATA_WRONLY(SELF%IODP)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IODP)
#else
         !$acc enter data attach(SELF%IODP)
#endif
         CALL SELF%F_IOBND%GET_DEVICE_DATA_WRONLY(SELF%IOBND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IOBND)
#else
         !$acc enter data attach(SELF%IOBND)
#endif
         CALL SELF%F_IBRMEM%GET_DEVICE_DATA_WRONLY(SELF%IBRMEM)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IBRMEM)
#else
         !$acc enter data attach(SELF%IBRMEM)
#endif
         CALL SELF%F_DELLAM1%GET_DEVICE_DATA_WRONLY(SELF%DELLAM1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%DELLAM1)
#else
         !$acc enter data attach(SELF%DELLAM1)
#endif
         CALL SELF%F_COSPHM1%GET_DEVICE_DATA_WRONLY(SELF%COSPHM1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%COSPHM1)
#else
         !$acc enter data attach(SELF%COSPHM1)
#endif
         CALL SELF%F_DEPTH%GET_DEVICE_DATA_WRONLY(SELF%DEPTH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%DEPTH)
#else
         !$acc enter data attach(SELF%DEPTH)
#endif
         CALL SELF%F_EMAXDPT%GET_DEVICE_DATA_WRONLY(SELF%EMAXDPT)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%EMAXDPT)
#else
         !$acc enter data attach(SELF%EMAXDPT)
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
      ENDIF

   END SUBROUTINE ENVIRONMENT_GET_DEVICE_DATA_WRONLY

   SUBROUTINE ENVIRONMENT_GET_DEVICE_DATA_RDONLY(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
             CALL SELF%F_INDEP%GET_DEVICE_DATA_RDONLY(SELF%INDEP)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%INDEP)
#else
             !$acc enter data attach(SELF%INDEP)
#endif
         ENDIF
         IF(L_IODP)THEN
             CALL SELF%F_IODP%GET_DEVICE_DATA_RDONLY(SELF%IODP)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IODP)
#else
             !$acc enter data attach(SELF%IODP)
#endif
         ENDIF
         IF(L_IOBND)THEN
             CALL SELF%F_IOBND%GET_DEVICE_DATA_RDONLY(SELF%IOBND)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IOBND)
#else
             !$acc enter data attach(SELF%IOBND)
#endif
         ENDIF
         IF(L_IBRMEM)THEN
             CALL SELF%F_IBRMEM%GET_DEVICE_DATA_RDONLY(SELF%IBRMEM)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%IBRMEM)
#else
             !$acc enter data attach(SELF%IBRMEM)
#endif
         ENDIF
         IF(L_DELLAM1)THEN
             CALL SELF%F_DELLAM1%GET_DEVICE_DATA_RDONLY(SELF%DELLAM1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%DELLAM1)
#else
             !$acc enter data attach(SELF%DELLAM1)
#endif
         ENDIF
         IF(L_COSPHM1)THEN
             CALL SELF%F_COSPHM1%GET_DEVICE_DATA_RDONLY(SELF%COSPHM1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%COSPHM1)
#else
             !$acc enter data attach(SELF%COSPHM1)
#endif
         ENDIF
         IF(L_DEPTH)THEN
             CALL SELF%F_DEPTH%GET_DEVICE_DATA_RDONLY(SELF%DEPTH)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%DEPTH)
#else
             !$acc enter data attach(SELF%DEPTH)
#endif
         ENDIF
         IF(L_EMAXDPT)THEN
             CALL SELF%F_EMAXDPT%GET_DEVICE_DATA_RDONLY(SELF%EMAXDPT)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%EMAXDPT)
#else
             !$acc enter data attach(SELF%EMAXDPT)
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
      ELSE
         !... copy entire struct
         CALL SELF%F_INDEP%GET_DEVICE_DATA_RDONLY(SELF%INDEP)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%INDEP)
#else
         !$acc enter data attach(SELF%INDEP)
#endif
         CALL SELF%F_IODP%GET_DEVICE_DATA_RDONLY(SELF%IODP)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IODP)
#else
         !$acc enter data attach(SELF%IODP)
#endif
         CALL SELF%F_IOBND%GET_DEVICE_DATA_RDONLY(SELF%IOBND)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IOBND)
#else
         !$acc enter data attach(SELF%IOBND)
#endif
         CALL SELF%F_IBRMEM%GET_DEVICE_DATA_RDONLY(SELF%IBRMEM)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%IBRMEM)
#else
         !$acc enter data attach(SELF%IBRMEM)
#endif
         CALL SELF%F_DELLAM1%GET_DEVICE_DATA_RDONLY(SELF%DELLAM1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%DELLAM1)
#else
         !$acc enter data attach(SELF%DELLAM1)
#endif
         CALL SELF%F_COSPHM1%GET_DEVICE_DATA_RDONLY(SELF%COSPHM1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%COSPHM1)
#else
         !$acc enter data attach(SELF%COSPHM1)
#endif
         CALL SELF%F_DEPTH%GET_DEVICE_DATA_RDONLY(SELF%DEPTH)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%DEPTH)
#else
         !$acc enter data attach(SELF%DEPTH)
#endif
         CALL SELF%F_EMAXDPT%GET_DEVICE_DATA_RDONLY(SELF%EMAXDPT)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%EMAXDPT)
#else
         !$acc enter data attach(SELF%EMAXDPT)
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
      ENDIF

   END SUBROUTINE ENVIRONMENT_GET_DEVICE_DATA_RDONLY

   SUBROUTINE ENVIRONMENT_SYNC_DEVICE_RDWR(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR, QUEUE)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
             CALL SELF%F_INDEP%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_IODP)THEN
             CALL SELF%F_IODP%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_IOBND)THEN
             CALL SELF%F_IOBND%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_IBRMEM)THEN
             CALL SELF%F_IBRMEM%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_DELLAM1)THEN
             CALL SELF%F_DELLAM1%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_COSPHM1)THEN
             CALL SELF%F_COSPHM1%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_DEPTH)THEN
             CALL SELF%F_DEPTH%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_EMAXDPT)THEN
             CALL SELF%F_EMAXDPT%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_INDEP%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_IODP%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_IOBND%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_IBRMEM%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_DELLAM1%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_COSPHM1%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_DEPTH%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_EMAXDPT%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE ENVIRONMENT_SYNC_DEVICE_RDWR

   SUBROUTINE ENVIRONMENT_SYNC_DEVICE_RDONLY(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR, QUEUE)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
             CALL SELF%F_INDEP%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_IODP)THEN
             CALL SELF%F_IODP%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_IOBND)THEN
             CALL SELF%F_IOBND%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_IBRMEM)THEN
             CALL SELF%F_IBRMEM%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_DELLAM1)THEN
             CALL SELF%F_DELLAM1%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_COSPHM1)THEN
             CALL SELF%F_COSPHM1%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_DEPTH)THEN
             CALL SELF%F_DEPTH%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_EMAXDPT)THEN
             CALL SELF%F_EMAXDPT%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_INDEP%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_IODP%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_IOBND%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_IBRMEM%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_DELLAM1%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_COSPHM1%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_DEPTH%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_EMAXDPT%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE ENVIRONMENT_SYNC_DEVICE_RDONLY

   SUBROUTINE ENVIRONMENT_GET_HOST_DATA_RDONLY(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%INDEP)
#else
             !$acc exit data detach(SELF%INDEP)
#endif
             CALL SELF%F_INDEP%GET_HOST_DATA_RDONLY(SELF%INDEP)
         ENDIF
         IF(L_IODP)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%IODP)
#else
             !$acc exit data detach(SELF%IODP)
#endif
             CALL SELF%F_IODP%GET_HOST_DATA_RDONLY(SELF%IODP)
         ENDIF
         IF(L_IOBND)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%IOBND)
#else
             !$acc exit data detach(SELF%IOBND)
#endif
             CALL SELF%F_IOBND%GET_HOST_DATA_RDONLY(SELF%IOBND)
         ENDIF
         IF(L_IBRMEM)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%IBRMEM)
#else
             !$acc exit data detach(SELF%IBRMEM)
#endif
             CALL SELF%F_IBRMEM%GET_HOST_DATA_RDONLY(SELF%IBRMEM)
         ENDIF
         IF(L_DELLAM1)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%DELLAM1)
#else
             !$acc exit data detach(SELF%DELLAM1)
#endif
             CALL SELF%F_DELLAM1%GET_HOST_DATA_RDONLY(SELF%DELLAM1)
         ENDIF
         IF(L_COSPHM1)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%COSPHM1)
#else
             !$acc exit data detach(SELF%COSPHM1)
#endif
             CALL SELF%F_COSPHM1%GET_HOST_DATA_RDONLY(SELF%COSPHM1)
         ENDIF
         IF(L_DEPTH)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%DEPTH)
#else
             !$acc exit data detach(SELF%DEPTH)
#endif
             CALL SELF%F_DEPTH%GET_HOST_DATA_RDONLY(SELF%DEPTH)
         ENDIF
         IF(L_EMAXDPT)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%EMAXDPT)
#else
             !$acc exit data detach(SELF%EMAXDPT)
#endif
             CALL SELF%F_EMAXDPT%GET_HOST_DATA_RDONLY(SELF%EMAXDPT)
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
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%INDEP)
#else
         !$acc exit data detach(SELF%INDEP)
#endif
         CALL SELF%F_INDEP%GET_HOST_DATA_RDONLY(SELF%INDEP)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%IODP)
#else
         !$acc exit data detach(SELF%IODP)
#endif
         CALL SELF%F_IODP%GET_HOST_DATA_RDONLY(SELF%IODP)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%IOBND)
#else
         !$acc exit data detach(SELF%IOBND)
#endif
         CALL SELF%F_IOBND%GET_HOST_DATA_RDONLY(SELF%IOBND)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%IBRMEM)
#else
         !$acc exit data detach(SELF%IBRMEM)
#endif
         CALL SELF%F_IBRMEM%GET_HOST_DATA_RDONLY(SELF%IBRMEM)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%DELLAM1)
#else
         !$acc exit data detach(SELF%DELLAM1)
#endif
         CALL SELF%F_DELLAM1%GET_HOST_DATA_RDONLY(SELF%DELLAM1)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%COSPHM1)
#else
         !$acc exit data detach(SELF%COSPHM1)
#endif
         CALL SELF%F_COSPHM1%GET_HOST_DATA_RDONLY(SELF%COSPHM1)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%DEPTH)
#else
         !$acc exit data detach(SELF%DEPTH)
#endif
         CALL SELF%F_DEPTH%GET_HOST_DATA_RDONLY(SELF%DEPTH)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%EMAXDPT)
#else
         !$acc exit data detach(SELF%EMAXDPT)
#endif
         CALL SELF%F_EMAXDPT%GET_HOST_DATA_RDONLY(SELF%EMAXDPT)
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
      ENDIF

   END SUBROUTINE ENVIRONMENT_GET_HOST_DATA_RDONLY

   SUBROUTINE ENVIRONMENT_GET_HOST_DATA_RDWR(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR

      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%INDEP)
#else
             !$acc exit data detach(SELF%INDEP)
#endif
             CALL SELF%F_INDEP%GET_HOST_DATA_RDWR(SELF%INDEP)
         ENDIF
         IF(L_IODP)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%IODP)
#else
             !$acc exit data detach(SELF%IODP)
#endif
             CALL SELF%F_IODP%GET_HOST_DATA_RDWR(SELF%IODP)
         ENDIF
         IF(L_IOBND)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%IOBND)
#else
             !$acc exit data detach(SELF%IOBND)
#endif
             CALL SELF%F_IOBND%GET_HOST_DATA_RDWR(SELF%IOBND)
         ENDIF
         IF(L_IBRMEM)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%IBRMEM)
#else
             !$acc exit data detach(SELF%IBRMEM)
#endif
             CALL SELF%F_IBRMEM%GET_HOST_DATA_RDWR(SELF%IBRMEM)
         ENDIF
         IF(L_DELLAM1)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%DELLAM1)
#else
             !$acc exit data detach(SELF%DELLAM1)
#endif
             CALL SELF%F_DELLAM1%GET_HOST_DATA_RDWR(SELF%DELLAM1)
         ENDIF
         IF(L_COSPHM1)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%COSPHM1)
#else
             !$acc exit data detach(SELF%COSPHM1)
#endif
             CALL SELF%F_COSPHM1%GET_HOST_DATA_RDWR(SELF%COSPHM1)
         ENDIF
         IF(L_DEPTH)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%DEPTH)
#else
             !$acc exit data detach(SELF%DEPTH)
#endif
             CALL SELF%F_DEPTH%GET_HOST_DATA_RDWR(SELF%DEPTH)
         ENDIF
         IF(L_EMAXDPT)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%EMAXDPT)
#else
             !$acc exit data detach(SELF%EMAXDPT)
#endif
             CALL SELF%F_EMAXDPT%GET_HOST_DATA_RDWR(SELF%EMAXDPT)
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
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%INDEP)
#else
         !$acc exit data detach(SELF%INDEP)
#endif
         CALL SELF%F_INDEP%GET_HOST_DATA_RDWR(SELF%INDEP)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%IODP)
#else
         !$acc exit data detach(SELF%IODP)
#endif
         CALL SELF%F_IODP%GET_HOST_DATA_RDWR(SELF%IODP)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%IOBND)
#else
         !$acc exit data detach(SELF%IOBND)
#endif
         CALL SELF%F_IOBND%GET_HOST_DATA_RDWR(SELF%IOBND)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%IBRMEM)
#else
         !$acc exit data detach(SELF%IBRMEM)
#endif
         CALL SELF%F_IBRMEM%GET_HOST_DATA_RDWR(SELF%IBRMEM)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%DELLAM1)
#else
         !$acc exit data detach(SELF%DELLAM1)
#endif
         CALL SELF%F_DELLAM1%GET_HOST_DATA_RDWR(SELF%DELLAM1)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%COSPHM1)
#else
         !$acc exit data detach(SELF%COSPHM1)
#endif
         CALL SELF%F_COSPHM1%GET_HOST_DATA_RDWR(SELF%COSPHM1)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%DEPTH)
#else
         !$acc exit data detach(SELF%DEPTH)
#endif
         CALL SELF%F_DEPTH%GET_HOST_DATA_RDWR(SELF%DEPTH)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%EMAXDPT)
#else
         !$acc exit data detach(SELF%EMAXDPT)
#endif
         CALL SELF%F_EMAXDPT%GET_HOST_DATA_RDWR(SELF%EMAXDPT)
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
      ENDIF

   END SUBROUTINE ENVIRONMENT_GET_HOST_DATA_RDWR

   SUBROUTINE ENVIRONMENT_SYNC_HOST_RDONLY(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR, QUEUE)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
             CALL SELF%F_INDEP%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_IODP)THEN
             CALL SELF%F_IODP%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_IOBND)THEN
             CALL SELF%F_IOBND%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_IBRMEM)THEN
             CALL SELF%F_IBRMEM%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_DELLAM1)THEN
             CALL SELF%F_DELLAM1%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_COSPHM1)THEN
             CALL SELF%F_COSPHM1%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_DEPTH)THEN
             CALL SELF%F_DEPTH%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_EMAXDPT)THEN
             CALL SELF%F_EMAXDPT%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_INDEP%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_IODP%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_IOBND%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_IBRMEM%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_DELLAM1%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_COSPHM1%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_DEPTH%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_EMAXDPT%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_HOST_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE ENVIRONMENT_SYNC_HOST_RDONLY

   SUBROUTINE ENVIRONMENT_SYNC_HOST_RDWR(SELF, INDEP, IODP, IOBND, IBRMEM, DELLAM1, COSPHM1, DEPTH, EMAXDPT, UCUR, VCUR, QUEUE)
      TYPE(ENVIRONMENT) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: INDEP
      LOGICAL :: L_INDEP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IODP
      LOGICAL :: L_IODP = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IOBND
      LOGICAL :: L_IOBND = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: IBRMEM
      LOGICAL :: L_IBRMEM = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DELLAM1
      LOGICAL :: L_DELLAM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: COSPHM1
      LOGICAL :: L_COSPHM1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: DEPTH
      LOGICAL :: L_DEPTH = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: EMAXDPT
      LOGICAL :: L_EMAXDPT = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: UCUR
      LOGICAL :: L_UCUR = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: VCUR
      LOGICAL :: L_VCUR = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(INDEP)) L_INDEP = INDEP
     IF(PRESENT(IODP)) L_IODP = IODP
     IF(PRESENT(IOBND)) L_IOBND = IOBND
     IF(PRESENT(IBRMEM)) L_IBRMEM = IBRMEM
     IF(PRESENT(DELLAM1)) L_DELLAM1 = DELLAM1
     IF(PRESENT(COSPHM1)) L_COSPHM1 = COSPHM1
     IF(PRESENT(DEPTH)) L_DEPTH = DEPTH
     IF(PRESENT(EMAXDPT)) L_EMAXDPT = EMAXDPT
     IF(PRESENT(UCUR)) L_UCUR = UCUR
     IF(PRESENT(VCUR)) L_VCUR = VCUR


      IF(L_INDEP.OR.L_IODP.OR.L_IOBND.OR.L_IBRMEM.OR.L_DELLAM1.OR.L_COSPHM1.OR.L_DEPTH.OR.L_EMAXDPT.OR.L_UCUR.OR.L_VCUR)THEN
         !... copy only selected members
         IF(L_INDEP)THEN
             CALL SELF%F_INDEP%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_IODP)THEN
             CALL SELF%F_IODP%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_IOBND)THEN
             CALL SELF%F_IOBND%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_IBRMEM)THEN
             CALL SELF%F_IBRMEM%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_DELLAM1)THEN
             CALL SELF%F_DELLAM1%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_COSPHM1)THEN
             CALL SELF%F_COSPHM1%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_DEPTH)THEN
             CALL SELF%F_DEPTH%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_EMAXDPT)THEN
             CALL SELF%F_EMAXDPT%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_UCUR)THEN
             CALL SELF%F_UCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_VCUR)THEN
             CALL SELF%F_VCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_INDEP%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_IODP%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_IOBND%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_IBRMEM%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_DELLAM1%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_COSPHM1%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_DEPTH%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_EMAXDPT%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_UCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_VCUR%SYNC_HOST_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE ENVIRONMENT_SYNC_HOST_RDWR

   SUBROUTINE ENVIRONMENT_DELETE_DEVICE_DATA(SELF)
      TYPE(ENVIRONMENT) :: SELF

      CALL SELF%F_INDEP%DELETE_DEVICE_DATA()
      CALL SELF%F_IODP%DELETE_DEVICE_DATA()
      CALL SELF%F_IOBND%DELETE_DEVICE_DATA()
      CALL SELF%F_IBRMEM%DELETE_DEVICE_DATA()
      CALL SELF%F_DELLAM1%DELETE_DEVICE_DATA()
      CALL SELF%F_COSPHM1%DELETE_DEVICE_DATA()
      CALL SELF%F_DEPTH%DELETE_DEVICE_DATA()
      CALL SELF%F_EMAXDPT%DELETE_DEVICE_DATA()
      CALL SELF%F_UCUR%DELETE_DEVICE_DATA()
      CALL SELF%F_VCUR%DELETE_DEVICE_DATA()

   END SUBROUTINE ENVIRONMENT_DELETE_DEVICE_DATA

END MODULE ENVIRONMENT_TYPE_MOD
