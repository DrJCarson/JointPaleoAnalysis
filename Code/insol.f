       ! Fortran subroutine for reading precession and obliquity
       ! Precession and obliquity are computed after Berger (1978)
       ! Precession is reconstructed as a d'Alambert series. This
       ! is not quite the proceduce advised by Berger, but is the most
       ! efficient for the purpose of investigating the dynamical
       ! response of a simple model

       ! NAP : Number of terms in the d'Alember decomp. of precession : 
       !       chose at least 20
       ! NAO : Same but for obliquity : at least 20
       
       ! USage : execute first 'READ_PRECESSION' and 'READ_OBLIQUITY'
       ! with the desired NAO and NAP. Theses routines output
       ! AMPlitudes, OMEgas (angular velocity) and phase (ANG)
       ! as arrays of dimension NAP and NAO, respectively

       ! SCALE TIME defaults to 10^4 : this means that time units
       ! are 10,000 years

       ! Then calles ASTRO, which returns a normalised version
       ! of PRECESSIOS (XPRE) and its derivative (DXPREDT) and
       ! obliquity (XOBL) with its derivate (DXOBLDT)

       ! Compile in FORTRAN90, e.g. : gfortan -c -o insol.o insol.f 
       ! or for use with R routine : R CMD SHLIB insol.f
  
       ! supplied : obliquity and precession must be in a 
       !            a Data directory.


       ! With this procedure a good approximation of 
       ! Summer Solstice Insolation at 65 degrees North
       ! may be obtained as follows (normalised)

       ! 0.7639 * XPRE + 0.4756 * XOBL       

       ! Author : M. Crucifix, Summer 2011. 


       SUBROUTINE READ_PRECESSION(NAP,AMPPRE,OMEPRE, ANGPRE)

       IMPLICIT NONE
       INTEGER NAP
       INTEGER I    ! COUNTER
       DOUBLE PRECISION AMPPRE(NAP),OMEPRE(NAP), ANGPRE(NAP)
       DOUBLE PRECISION TMP1,TMP2,TMP3 ! TMP
       DOUBLE PRECISION SCALEPREC,SCALETIME
       DOUBLE PRECISION TWOPI, TWOPI360
       PARAMETER (SCALETIME = 1D4 )
       PARAMETER (SCALEPREC = 1./0.01864561223205 )
       PARAMETER (TWOPI     = 2*3.1415926535897   )
       PARAMETER (TWOPI360  = TWOPI / 360.        )

       OPEN (10, FILE='Data/precession.dat'   ,STATUS='OLD')
       DO I=1,NAP
         READ (10,*) TMP1,TMP2,TMP3
         AMPPRE(I)=TMP1*SCALEPREC
         OMEPRE(I)=TMP2*(TWOPI360/(60.*60.))*SCALETIME
         ANGPRE(I)=TMP3*TWOPI360
       ENDDO
       CLOSE(10)
       END

      SUBROUTINE READ_OBLIQUITY(NAO,AMPOBL,OMEOBL, ANGOBL)

       IMPLICIT NONE
       INTEGER NAO
       INTEGER I    ! COUNTER
       DOUBLE PRECISION AMPOBL(NAO),OMEOBL(NAO), ANGOBL(NAO)
       DOUBLE PRECISION TMP1,TMP2,TMP3 ! TMP
       DOUBLE PRECISION SCALEOBL,SCALETIME
       DOUBLE PRECISION TWOPI, TWOPI360
       PARAMETER (SCALETIME = 1D4 )
       PARAMETER (SCALEOBL = 1./2462.2214466 )
       PARAMETER (TWOPI     = 2*3.1415926535897   )
       PARAMETER (TWOPI360  = TWOPI / 360.        )

      OPEN (10, FILE='Data/obliquity.dat'   ,STATUS='OLD')
       DO I=1,NAO
         READ (10,*) TMP1,TMP2,TMP3
         AMPOBL(I)=TMP1*SCALEOBL
         OMEOBL(I)=TMP2*(TWOPI360/(60.*60.))*SCALETIME
         ANGOBL(I)=TMP3*TWOPI360
       ENDDO
       CLOSE(10)
       END


       SUBROUTINE ASTRO(NAP, NAO, T,AMPPRE, OMEPRE, ANGPRE, AMPOBL, 
     c                OMEOBL, ANGOBL, XPRE, DXPREDT, XOBL, DXOBLDT)
        IMPLICIT NONE
        INTEGER  NAP, NAO
        DOUBLE PRECISION T,AMPPRE(NAP), OMEPRE(NAP), ANGPRE(NAP)  ! INPUT
        DOUBLE PRECISION   AMPOBL(NAO), OMEOBL(NAO), ANGOBL(NAO)  ! INPUT
        DOUBLE PRECISION XPRE, DXPREDT                 ! OUTPUT
        DOUBLE PRECISION XOBL, DXOBLDT                 ! OUTPUT
        DOUBLE PRECISION SOMEP(NAP),COMEP(NAP)             ! TMP
        DOUBLE PRECISION SOMEO(NAO),COMEO(NAO)             ! TMP

        ! PRECESSION
        SOMEP = DSIN(OMEPRE*T + ANGPRE)
        COMEP = DCOS(OMEPRE*T + ANGPRE)
        XPRE  = SUM(AMPPRE*SOMEP)
        DXPREDT = SUM(OMEPRE*(AMPPRE*COMEP))

        ! OBLIQUITY
        SOMEO = DSIN(OMEOBL*T + ANGOBL)
        COMEO = DCOS(OMEOBL*T + ANGOBL)
        XOBL  = SUM(AMPOBL*COMEO)
        DXOBLDT = - SUM(OMEOBL*(AMPOBL*SOMEO))
      END 


