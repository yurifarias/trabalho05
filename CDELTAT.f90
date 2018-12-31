SUBROUTINE CDELTAT(NELE, NNOS, CELE, CNOS, LAR, ALT, MEL, MESP, DT)
    ! Calcular DELTAT para problema transiente
    IMPLICIT NONE

    INTEGER, INTENT(in) :: NELE, NNOS
    INTEGER, DIMENSION(NELE,2), INTENT(in) :: CELE
    REAL(kind=8), DIMENSION(NNOS,2), INTENT(in) :: CNOS
    REAL(kind=8), DIMENSION(NELE), INTENT(in) :: LAR, ALT, MEL, MESP
    REAL(kind=8), INTENT(out) :: DT

    INTEGER :: i
    REAL(kind=16), PARAMETER :: PI_16 = 4 * ATAN(1.0_16)
    REAL(kind=8), DIMENSION(NELE) :: area, iner, comp, freq_nat
    REAL(kind=8), DIMENSION(NELE, 2) :: proj

    DO i=1, NELE
        IF (LAR(i) == 0) THEN
            area(i) = PI_16 * (ALT(i) ** 2) / 4
            iner(i) = PI_16 * (ALT(i) ** 4) / 64
        ELSE
            area(i) = LAR(i) * ALT(i)
            iner(i) = LAR(i) * ALT(i) ** 3 / 12
        END IF
        proj(i,1) = CNOS(CELE(i,2), 1) - CNOS(CELE(i,1), 1)         ! Projeção sobre o eixo x
        proj(i,2) = CNOS(CELE(i,2), 2) - CNOS(CELE(i,1), 2)         ! Projeção sobre o eixo y
        comp(i) = ((proj(i,1) ** 2) + (proj(i,2) ** 2)) ** 0.5      ! Comprimento
    END DO

    DO i=1, NELE
        freq_nat(i) = (((PI_16 ** 2) / (comp(i) ** 2))) * ((MEL(i) * iner(i)) / (MESP(i) * area(i))) ** 0.5
    END DO

    DT = 2 / (MAXVAL(freq_nat) ** 2)

END SUBROUTINE CDELTAT