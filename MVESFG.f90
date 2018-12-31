SUBROUTINE MVESFG(NELE, NNOS, CELE, ROTS, ESFL, CARN, ESFG)
    ! Subrotina MVESFL: Montar Vetor de ESForços Globais da estrutura
    IMPLICIT NONE

    INTEGER, INTENT(in) :: NELE, NNOS
    INTEGER, DIMENSION(NELE,2), INTENT(in) :: CELE
    REAL(kind=8), DIMENSION(6,6,NELE), INTENT(in) :: ROTS
    REAL(kind=8), DIMENSION(6,NELE), INTENT(in) :: ESFL
    REAL(kind=8), DIMENSION(3*NNOS), INTENT(in) :: CARN
    REAL(kind=8), DIMENSION(3*NNOS), INTENT(out) :: ESFG

    INTEGER :: i, j, k, noj, nok
    REAL(kind=8), DIMENSION(6) :: vesf_loc, vesf_glo
    REAL(kind=8), DIMENSION(6,6) :: mrot

    ! NELE: variável que armazena número de elementos
    ! NNOS: variável que armazena número de nós
    ! CELE: matriz que armazena índices dos nós iniciais e finais de cada elemento [nói,nóf]
    ! ROTS: matriz 3D que armazena a matriz de rotação local de cada elemento
    ! ESFL: matriz 2D que armazena vetor de esforço local em cada elemento [6,num_elem]
    ! CARN: vetor que armazena esforços nodais
    ! ESFG: vetor que armazena esforços globais na estrutura

    ! i: variável para ser iterada
    ! j: variável para ser iterada
    ! k: variável para ser iterada
    ! noj: variável de auxílio para montagem da matriz de rigidez da estrutura (define linha)
    ! nok: variável de auxílio para montagem da matriz de rigidez da estrutura (define coluna)
    ! vesf_loc: vetor auxiliar para armazenar vetor de esforços locais de cada elemento
    ! vesf_glo: vetor auxiliar para armazenar vetor de esforços globais de cada elemento
    ! mrot: matriz auxiliar para armazenar matriz de rotação de cada elemento

    ESFG = 0.0

    DO i=1, NELE
        ! Montagem do vetor de esforços locais e matriz de rotação do elemento i
        DO j=1, 6
            vesf_loc(j) = ESFL(j,i)                 ! Vetor de esforços locais
            DO k=1, 6
                mrot(j,k) = ROTS(j,k,i)             ! Matriz de rotação
            END DO
        END DO

        ! Cálculo do vetor de esforços globais {esfg} = [ROT]t{esfl}
        vesf_glo = MATMUL(TRANSPOSE(mrot),vesf_loc)

        ! Endereçamento do vetor de solicitações globais para vetor esforços globais
        DO j=-2, 0
            noj = 3 * CELE(i,1) + j
            nok = 3 * CELE(i,2) + j
            ESFG(noj) = ESFG(noj) + vesf_glo(3 + j)
            ESFG(nok) = ESFG(nok) + vesf_glo(6 + j)
        ENDDO
    END DO

    ESFG = ESFG + CARN

END SUBROUTINE MVESFG