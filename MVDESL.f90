SUBROUTINE MVDESL(NELE, NNOS, CELE, DESG, ROTS, DESL)
    ! Subrotina MVDESL: Montar Vetor de DESlocamentos Locais de cada elemento
    IMPLICiT NONE

    INTEGER, INTENT(in) :: NELE, NNOS
    INTEGER, DIMENSION(NELE,2), INTENT(in) :: CELE
    REAL(kind=8), DIMENSION(3*NNOS), INTENT(in) :: DESG
    REAL(kind=8), DIMENSION(6,6,NELE), INTENT(in) :: ROTS
    REAL(kind=8), DIMENSION(6,NELE), INTENT(out) :: DESL

    INTEGER :: i, j, k, noj, nok
    REAL(kind=8), DIMENSION(6) :: vdes_glo, vdes_loc
    REAL(kind=8), DIMENSION(6,6) :: mrot

    ! NELE: variável que armazena número de elementos
    ! NNOS: variável que armazena número de nós
    ! DESG: vetor que armazena deslocamentos globais na estrutura
    ! CELE: matriz que armazena índices dos nós iniciais e finais de cada elemento [nói,nóf]
    ! ROTS: matriz 3D que armazena a matriz de rotação local de cada elemento
    ! DESL: matriz 2D que armazena vetor de deslocamentos locais em cada elemento [6,num_elem]

    ! i: variável para ser iterada
    ! j: variável para ser iterada
    ! noj: variável de auxílio para montagem da matriz de rigidez da estrutura (define linha)
    ! nok: variável de auxílio para montagem da matriz de rigidez da estrutura (define coluna)
    ! vdes_glo: vetor auxiliar para armazenar vetor de deslocamentos globais de cada elemento
    ! vdes_loc: vetor auxiliar para armazenar vetor de deslocamentos locais de cada elemento
    ! mrot: matriz auxiliar para armazenar matriz de rotação de cada elemento

    DO i=1, NELE
        ! Endereçamento do vetor de deslocamentos globais da estrutura para vetor
        ! de deslocamentos globais do elemento
        DO j=-2, 0
            noj = 3 * CELE(i,1) + j             ! Linha referente ao nó inicial
            nok = 3 * CELE(i,2) + j             ! Linha referente ao nó final
            vdes_glo(3+j) = DESG(noj)
            vdes_glo(6+j) = DESG(nok)
        ENDDO

        ! Montagem do vetor de deslocamentos globais e matriz de rotação do elemento i
        DO j=1, 6
            DO k=1, 6
                mrot(j,k) = ROTS(j,k,i)         ! Matriz de rotação
            END DO
        END DO

        ! Cálculo do vetor de esforços globais {desl} = [ROT]{desg}
        vdes_loc = MATMUL(mrot,vdes_glo)

        ! Montar matriz 2D para armazenar deslocamentos locais de cada elemento
        DO j=1, 6
            DESL(j,i) = vdes_loc(j)
        END DO
    END DO
END SUBROUTINE MVDESL