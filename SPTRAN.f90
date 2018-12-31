SUBROUTINE SPTRAN(NNOS, RIGG, MASG, AMORT, CARM, DT)
    ! Subrotina SPTRAN: Solucionar Problema TRANsiente
    IMPLICIT NONE

    INTEGER, INTENT(in) :: NNOS
    REAL(kind=8), INTENT(in) :: DT
    REAL(kind=8), DIMENSION(3*NNOS) :: CARM
    REAL(kind=8), DIMENSION(3*NNOS,3*NNOS), INTENT(in) :: RIGG, MASG, AMORT

    INTEGER :: n, i, tf = 60
    INTEGER, PARAMETER :: f10 = 110
    REAL(kind=8) :: ti
    REAL(kind=8), DIMENSION(3*NNOS) :: U_zero_dp, U_menos_um, Un_menos_um, Un, Un_mais_um, Fchap_n_mais_um
    REAL(kind=8), DIMENSION(3*NNOS) :: U_zero, U_zero_p
    REAL(kind=8), DIMENSION(3*NNOS,3*NNOS) :: R2, R3, R4, aux1, aux2

    ! NNOS: variável que armazena número de nós
    ! RIGG: matriz de rigidez global da estrutura
    ! MASG: matriz de massa global da estrutura
    ! AMORT: matriz de amortecimento da estrutura
    ! CARM: vetor de cargas móveis nodais

    ! delta_t: espaço de tempo para incremento
    ! U_zero: vetor deslocamento no tempo inicial (ti = 0)
    ! U_zero_p: vetor velocidade no tempo inicial (ti = 0)
    ! U_zero_dp: vetor aceleração no tempo inicial (ti = 0)
    ! U_menos_um: vetor deslocamento no tempo t-1 (t-1 = ti - delta_t)

    ! Cálculo das acelerações nodais iniciais U_zero_dp (no momento t = 0)
    CALL DLSLRGa(3*NNOS,3*NNOS,MASG,CARM,1,U_zero_dp)

    U_zero = 0
    U_zero_p = 0

    ! Cálculo do deslocamento U_menos_um
    U_menos_um = 0.5 * U_zero_dp * (DT ** 2)

    ! Cálculo de R2, R3 e R4
    aux1 = AMORT / (2 * DT)
    aux2 = MASG / (DT ** 2)
    R2 = aux1 + aux2                                        ! C/2*dt + M/dt**2
    R3 = aux2 - aux1                                        ! M/dt**2 - C/2*dt
    R4 = RIGG - 2 * aux2                                    ! K - 2*M/dt**2
    ti = 0
    n = 0

    Un_menos_um = U_menos_um
    Un = U_zero

    OPEN(f10,file='respostatransiente.csv')

    WRITE(f10,'(*(g0))') 'ITERACAO,TEMPO (s),DESLOCAMENTO (m)'
    WRITE(f10,'(*(g0))') n, ',',ti, ',',Un(14)

    DO WHILE (ti <= tf)
        IF (ti > 10) THEN
            CARM = 0
        END IF

        Fchap_n_mais_um = CARM - MATMUL(R4,Un) - MATMUL(R3,Un_menos_um)

        CALL DLSLRGa(3*NNOS,3*NNOS,R2,Fchap_n_mais_um,1,Un_mais_um)

        n = n + 1
        ti = n * DT

        IF (MOD(n,1000) == 0) THEN
            WRITE(f10,'(*(g0))') n, ',',ti, ',',Un_mais_um(14)
        END IF

        Un_menos_um = Un
        Un = Un_mais_um

    END DO

    CLOSE(f10)

END SUBROUTINE SPTRAN