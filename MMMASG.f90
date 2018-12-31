SUBROUTINE MMMASG(NELE, NNOS, CELE, MASL, ROTS, MASG)
    ! Subrotina MMMASG: Montar Matriz de MASsa Global da estrutura
    IMPLICIT NONE

    INTEGER, INTENT(in) :: NELE, NNOS
    INTEGER, DIMENSION(NELE,2) :: CELE
    REAL(kind=8), DIMENSION(6,6,NELE), INTENT(in) :: MASL, ROTS
    REAL(kind=8), DIMENSION(3*NNOS,3*NNOS), INTENT(out) :: MASG

    INTEGER :: i, j, k, noj, nok
    REAL(kind=8), DIMENSION(6,6) :: mmas_loc, mmas_glo, mrot

    ! NELE: variável que armazena número de elementos
    ! NNOS: variável que armazena número de nós
    ! CELE: matriz que armazena índices dos nós iniciais e finais de cada elemento [nói, nóf]
    ! MASL: matriz 3D que armazena matrizes de massa local de cada elemento
    ! ROTS: matriz 3D que armazena matrizes de rotação local de cada elemento
    ! MASG: matriz de massa global da estrutura

    ! i: variável para ser iterada
    ! j: variável para ser iterada
    ! k: variável para ser iterada
    ! noj: variável de auxílio para montagem da matriz de massa da estrutura (define linha)
    ! nok: variável de auxílio para montagem da matriz de massa da estrutura (define coluna)
    ! mmas_loc: matriz auxiliar para armazenar matriz de massa local de cada elemento
    ! mmas_glo: matriz auxiliar para armazenar matriz de massa global de cada elemento
    ! mrot: matriz auxiliar para armazenar matriz de rotação de cada elemento

    MASG = 0.0

    DO i=1, NELE
        ! Montagem da matriz de massa local e matriz de rotação do elemento i
        DO j=1, 6
            DO k=1, 6
                mmas_loc(j,k) = MASL(j,k,i)         ! Matriz de massa local
                mrot(j,k) = ROTS(j,k,i)             ! Matriz de rotação
            END DO
        END DO

        ! Cálculo da matriz de massa global do elemento i [Kg] = [R]t[Kl][R]
        mmas_glo = MATMUL(MATMUL(TRANSPOSE(mrot),mmas_loc),mrot)

        ! Início da montagem da matriz de massa global da estrutura

        ! Para (nó inicial, nó inicial)
        DO j=-2, 0
            noj = 3 * CELE(i,1) + j             ! Linha da matriz de massa global da estrutura
            DO k=-2, 0
                nok = 3 * CELE(i,1) + k         ! Coluna da matriz de massa global da estrutura
                MASG(noj, nok) = MASG(noj, nok) + mmas_glo(3 + j, 3 + k)
            ENDDO
        ENDDO

        ! Para (nó inicial, nó final)
        DO j=-2, 0
            noj = 3 * CELE(i,1) + j             ! Linha da matriz de massa global da estrutura
            DO k=-2, 0
                nok = 3 * CELE(i,2) + k         ! Coluna da matriz de massa global da estrutura
                MASG(noj, nok) = MASG(noj, nok) + mmas_glo(3 + j, 6 + k)
            ENDDO
        ENDDO

        ! Para (nó final, nó inicial)
        DO j=-2, 0
            noj = 3 * CELE(i,2) + j             ! Linha da matriz de massa global da estrutura
            DO k=-2, 0
                nok = 3 * CELE(i,1) + k         ! Coluna da matriz de massa global da estrutura
                MASG(noj, nok) = MASG(noj, nok) + mmas_glo(6 + j, 3 + k)
            ENDDO
        ENDDO

        ! Para (nó final, nó final)
        DO j=-2, 0
            noj = 3 * CELE(i,2) + j             ! Linha da matriz de massa global da estrutura
            DO k=-2, 0
                nok = 3 * CELE(i,2) + k         ! Coluna da matriz de massa global da estrutura
                MASG(noj, nok) = MASG(noj, nok) + mmas_glo(6 + j, 6 + k)
            ENDDO
        ENDDO

        ! Fim da montagem da matriz global de massa da estrutura
    END DO
END SUBROUTINE MMMASG