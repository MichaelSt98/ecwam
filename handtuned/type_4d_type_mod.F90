MODULE TYPE_4D_TYPE_MOD

   USE PARKIND_WAVE, ONLY : JWRB, JWIM, JWRO
   USE FIELD_MODULE, ONLY : FIELD_3RB, FIELD_2IM, FIELD_2RB, FIELD_3IM, FIELD_4RB, FIELD_4IM, FIELD_2RD, &
   &                        FIELD_2RM, FIELD_1IM, FIELD_1RB
   USE FIELD_FACTORY_MODULE, ONLY : FIELD_NEW, FIELD_DELETE
   IMPLICIT NONE

   PRIVATE

   INTERFACE ALLOC
      MODULE PROCEDURE TYPE_4D_ALLOC
   END INTERFACE ALLOC
   INTERFACE DEALLOC
      MODULE PROCEDURE TYPE_4D_DEALLOC
   END INTERFACE DEALLOC
   INTERFACE SYNC_DEVICE_RDWR
      MODULE PROCEDURE TYPE_4D_SYNC_DEVICE_RDWR
   END INTERFACE SYNC_DEVICE_RDWR
   INTERFACE SYNC_DEVICE_RDONLY
      MODULE PROCEDURE TYPE_4D_SYNC_DEVICE_RDONLY
   END INTERFACE SYNC_DEVICE_RDONLY
   INTERFACE SYNC_HOST_RDWR
      MODULE PROCEDURE TYPE_4D_SYNC_HOST_RDWR
   END INTERFACE SYNC_HOST_RDWR
   INTERFACE SYNC_HOST_RDONLY
      MODULE PROCEDURE TYPE_4D_SYNC_HOST_RDONLY
   END INTERFACE SYNC_HOST_RDONLY
   INTERFACE GET_DEVICE_DATA_RDWR
      MODULE PROCEDURE TYPE_4D_GET_DEVICE_DATA_RDWR
   END INTERFACE GET_DEVICE_DATA_RDWR
   INTERFACE GET_DEVICE_DATA_RDONLY
      MODULE PROCEDURE TYPE_4D_GET_DEVICE_DATA_RDONLY
   END INTERFACE GET_DEVICE_DATA_RDONLY
   INTERFACE GET_DEVICE_DATA_WRONLY
      MODULE PROCEDURE TYPE_4D_GET_DEVICE_DATA_WRONLY
   END INTERFACE GET_DEVICE_DATA_WRONLY
   INTERFACE GET_HOST_DATA_RDWR
      MODULE PROCEDURE TYPE_4D_GET_HOST_DATA_RDWR
   END INTERFACE GET_HOST_DATA_RDWR
   INTERFACE GET_HOST_DATA_RDONLY
      MODULE PROCEDURE TYPE_4D_GET_HOST_DATA_RDONLY
   END INTERFACE GET_HOST_DATA_RDONLY
   INTERFACE DELETE_DEVICE_DATA
      MODULE PROCEDURE TYPE_4D_DELETE_DEVICE_DATA
   END INTERFACE DELETE_DEVICE_DATA

   TYPE TYPE_4D
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_4RB), POINTER :: F_FL1 => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: FL1(:,:,:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_4RB), POINTER :: F_FL1 => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: FL1(:,:,:,:) => NULL()
#endif
#ifndef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_4RB), POINTER :: F_XLLWS => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XLLWS(:,:,:,:) => NULL()
#endif
#ifdef PARKIND1_SINGLE_NEMO
      CLASS(FIELD_4RB), POINTER :: F_XLLWS => NULL()
      REAL(KIND=JWRB), POINTER, CONTIGUOUS :: XLLWS(:,:,:,:) => NULL()
#endif
      LOGICAL :: LALLOC = .FALSE.
   END TYPE TYPE_4D

   PUBLIC :: TYPE_4D, ALLOC, DEALLOC
   PUBLIC :: SYNC_DEVICE_RDWR, SYNC_DEVICE_RDONLY, SYNC_HOST_RDWR, SYNC_HOST_RDONLY, GET_DEVICE_DATA_RDWR, GET_DEVICE_DATA_RDONLY,&
       & GET_DEVICE_DATA_WRONLY, GET_HOST_DATA_RDWR, GET_HOST_DATA_RDONLY, DELETE_DEVICE_DATA

   CONTAINS

   SUBROUTINE TYPE_4D_ALLOC(SELF, UBOUNDS, LBOUNDS)
      TYPE(TYPE_4D) :: SELF
      INTEGER(KIND=JWIM), INTENT(IN) :: UBOUNDS(4)
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: LBOUNDS(4)
      INTEGER(KIND=JWIM) :: LLBOUNDS(4)

      LLBOUNDS(:) = 1
      IF(PRESENT(LBOUNDS)) LLBOUNDS = LBOUNDS

      CALL FIELD_NEW(SELF%F_FL1, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_FL1%GET_HOST_DATA_RDWR(SELF%FL1)
      SELF%FL1(:,:,:,:) = 0_JWRB
      CALL FIELD_NEW(SELF%F_XLLWS, LBOUNDS=LLBOUNDS, UBOUNDS=UBOUNDS, PERSISTENT=.TRUE.)
      CALL SELF%F_XLLWS%GET_HOST_DATA_RDWR(SELF%XLLWS)
      SELF%XLLWS(:,:,:,:) = 0_JWRB

      SELF%LALLOC = .TRUE.
#ifdef OMPGPU
      !$omp target enter data map(to:SELF)
#else
      !$acc enter data copyin(SELF)
#endif

   END SUBROUTINE TYPE_4D_ALLOC

   SUBROUTINE TYPE_4D_DEALLOC(SELF)
      TYPE(TYPE_4D) :: SELF

      NULLIFY(SELF%FL1)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%FL1)
#else
      !$acc exit data detach(SELF%FL1)
#endif
      CALL FIELD_DELETE(SELF%F_FL1)
      NULLIFY(SELF%F_FL1)
      NULLIFY(SELF%XLLWS)
#ifdef OMPGPU
      !$omp target exit data map(release:SELF%XLLWS)
#else
      !$acc exit data detach(SELF%XLLWS)
#endif
      CALL FIELD_DELETE(SELF%F_XLLWS)
      NULLIFY(SELF%F_XLLWS)

      SELF%LALLOC = .FALSE.
#ifdef OMPGPU
      !$omp target exit data map(delete:SELF)
#else
      !$acc exit data delete(SELF)
#endif

   END SUBROUTINE TYPE_4D_DEALLOC

   SUBROUTINE TYPE_4D_GET_DEVICE_DATA_RDWR(SELF, FL1, XLLWS)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
             CALL SELF%F_FL1%GET_DEVICE_DATA_RDWR(SELF%FL1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%FL1)
#else
             !$acc enter data attach(SELF%FL1)
#endif
         ENDIF
         IF(L_XLLWS)THEN
             CALL SELF%F_XLLWS%GET_DEVICE_DATA_RDWR(SELF%XLLWS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%XLLWS)
#else
             !$acc enter data attach(SELF%XLLWS)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_FL1%GET_DEVICE_DATA_RDWR(SELF%FL1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%FL1)
#else
         !$acc enter data attach(SELF%FL1)
#endif
         CALL SELF%F_XLLWS%GET_DEVICE_DATA_RDWR(SELF%XLLWS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%XLLWS)
#else
         !$acc enter data attach(SELF%XLLWS)
#endif
      ENDIF

   END SUBROUTINE TYPE_4D_GET_DEVICE_DATA_RDWR

   SUBROUTINE TYPE_4D_GET_DEVICE_DATA_WRONLY(SELF, FL1, XLLWS)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
             CALL SELF%F_FL1%GET_DEVICE_DATA_WRONLY(SELF%FL1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%FL1)
#else
             !$acc enter data attach(SELF%FL1)
#endif
         ENDIF
         IF(L_XLLWS)THEN
             CALL SELF%F_XLLWS%GET_DEVICE_DATA_WRONLY(SELF%XLLWS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%XLLWS)
#else
             !$acc enter data attach(SELF%XLLWS)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_FL1%GET_DEVICE_DATA_WRONLY(SELF%FL1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%FL1)
#else
         !$acc enter data attach(SELF%FL1)
#endif
         CALL SELF%F_XLLWS%GET_DEVICE_DATA_WRONLY(SELF%XLLWS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%XLLWS)
#else
         !$acc enter data attach(SELF%XLLWS)
#endif
      ENDIF

   END SUBROUTINE TYPE_4D_GET_DEVICE_DATA_WRONLY

   SUBROUTINE TYPE_4D_GET_DEVICE_DATA_RDONLY(SELF, FL1, XLLWS)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
             CALL SELF%F_FL1%GET_DEVICE_DATA_RDONLY(SELF%FL1)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%FL1)
#else
             !$acc enter data attach(SELF%FL1)
#endif
         ENDIF
         IF(L_XLLWS)THEN
             CALL SELF%F_XLLWS%GET_DEVICE_DATA_RDONLY(SELF%XLLWS)
#ifdef OMPGPU
             !$omp target enter data map(to:SELF%XLLWS)
#else
             !$acc enter data attach(SELF%XLLWS)
#endif
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_FL1%GET_DEVICE_DATA_RDONLY(SELF%FL1)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%FL1)
#else
         !$acc enter data attach(SELF%FL1)
#endif
         CALL SELF%F_XLLWS%GET_DEVICE_DATA_RDONLY(SELF%XLLWS)
#ifdef OMPGPU
         !$omp target enter data map(to:SELF%XLLWS)
#else
         !$acc enter data attach(SELF%XLLWS)
#endif
      ENDIF

   END SUBROUTINE TYPE_4D_GET_DEVICE_DATA_RDONLY

   SUBROUTINE TYPE_4D_SYNC_DEVICE_RDWR(SELF, FL1, XLLWS, QUEUE)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
             CALL SELF%F_FL1%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_XLLWS)THEN
             CALL SELF%F_XLLWS%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_FL1%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
         CALL SELF%F_XLLWS%SYNC_DEVICE_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE TYPE_4D_SYNC_DEVICE_RDWR

   SUBROUTINE TYPE_4D_SYNC_DEVICE_RDONLY(SELF, FL1, XLLWS, QUEUE)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
             CALL SELF%F_FL1%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_XLLWS)THEN
             CALL SELF%F_XLLWS%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_FL1%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_XLLWS%SYNC_DEVICE_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE TYPE_4D_SYNC_DEVICE_RDONLY

   SUBROUTINE TYPE_4D_GET_HOST_DATA_RDONLY(SELF, FL1, XLLWS)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%FL1)
#else
             !$acc exit data detach(SELF%FL1)
#endif
             CALL SELF%F_FL1%GET_HOST_DATA_RDONLY(SELF%FL1)
         ENDIF
         IF(L_XLLWS)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%XLLWS)
#else
             !$acc exit data detach(SELF%XLLWS)
#endif
             CALL SELF%F_XLLWS%GET_HOST_DATA_RDONLY(SELF%XLLWS)
         ENDIF
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%FL1)
#else
         !$acc exit data detach(SELF%FL1)
#endif
         CALL SELF%F_FL1%GET_HOST_DATA_RDONLY(SELF%FL1)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%XLLWS)
#else
         !$acc exit data detach(SELF%XLLWS)
#endif
         CALL SELF%F_XLLWS%GET_HOST_DATA_RDONLY(SELF%XLLWS)
      ENDIF

   END SUBROUTINE TYPE_4D_GET_HOST_DATA_RDONLY

   SUBROUTINE TYPE_4D_GET_HOST_DATA_RDWR(SELF, FL1, XLLWS)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS

      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%FL1)
#else
             !$acc exit data detach(SELF%FL1)
#endif
             CALL SELF%F_FL1%GET_HOST_DATA_RDWR(SELF%FL1)
         ENDIF
         IF(L_XLLWS)THEN
#ifdef OMPGPU
             !$omp target exit data map(release:SELF%XLLWS)
#else
             !$acc exit data detach(SELF%XLLWS)
#endif
             CALL SELF%F_XLLWS%GET_HOST_DATA_RDWR(SELF%XLLWS)
         ENDIF
      ELSE
         !... copy entire struct
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%FL1)
#else
         !$acc exit data detach(SELF%FL1)
#endif
         CALL SELF%F_FL1%GET_HOST_DATA_RDWR(SELF%FL1)
#ifdef OMPGPU
         !$omp target exit data map(release:SELF%XLLWS)
#else
         !$acc exit data detach(SELF%XLLWS)
#endif
         CALL SELF%F_XLLWS%GET_HOST_DATA_RDWR(SELF%XLLWS)
      ENDIF

   END SUBROUTINE TYPE_4D_GET_HOST_DATA_RDWR

   SUBROUTINE TYPE_4D_SYNC_HOST_RDONLY(SELF, FL1, XLLWS, QUEUE)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
             CALL SELF%F_FL1%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
         IF(L_XLLWS)THEN
             CALL SELF%F_XLLWS%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_FL1%SYNC_HOST_RDONLY(QUEUE=QUEUE)
         CALL SELF%F_XLLWS%SYNC_HOST_RDONLY(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE TYPE_4D_SYNC_HOST_RDONLY

   SUBROUTINE TYPE_4D_SYNC_HOST_RDWR(SELF, FL1, XLLWS, QUEUE)
      TYPE(TYPE_4D) :: SELF
      LOGICAL, INTENT(IN), OPTIONAL :: FL1
      LOGICAL :: L_FL1 = .FALSE.
      LOGICAL, INTENT(IN), OPTIONAL :: XLLWS
      LOGICAL :: L_XLLWS = .FALSE.
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: QUEUE

     IF(PRESENT(FL1)) L_FL1 = FL1
     IF(PRESENT(XLLWS)) L_XLLWS = XLLWS


      IF(L_FL1.OR.L_XLLWS)THEN
         !... copy only selected members
         IF(L_FL1)THEN
             CALL SELF%F_FL1%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
         IF(L_XLLWS)THEN
             CALL SELF%F_XLLWS%SYNC_HOST_RDWR(QUEUE=QUEUE)
         ENDIF
      ELSE
         !... copy entire struct
         CALL SELF%F_FL1%SYNC_HOST_RDWR(QUEUE=QUEUE)
         CALL SELF%F_XLLWS%SYNC_HOST_RDWR(QUEUE=QUEUE)
      ENDIF

   END SUBROUTINE TYPE_4D_SYNC_HOST_RDWR

   SUBROUTINE TYPE_4D_DELETE_DEVICE_DATA(SELF)
      TYPE(TYPE_4D) :: SELF

      CALL SELF%F_FL1%DELETE_DEVICE_DATA()
      CALL SELF%F_XLLWS%DELETE_DEVICE_DATA()

   END SUBROUTINE TYPE_4D_DELETE_DEVICE_DATA

END MODULE TYPE_4D_TYPE_MOD
