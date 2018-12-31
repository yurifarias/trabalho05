SUBROUTINE MMROT(NELE, NNOS, CELE, CNOS, ROTS)
    ! Subrotina MMROT: Montar Matriz de ROTação de cada elemento
    IMPLICIT NONE

    INTEGER, INTENT(in) :: NELE, NNOS
    INTEGER, DIMENSION(NELE,2), INTENT(in) :: CELE
    REAL(kind=8), DIMENSION(NNOS,2), INTENT(in) :: CNOS
    REAL(kind=8), DIMENSION(6,6,NELE), INTENT(out) :: ROTS

    INTEGER :: i
    REAL(kind=8) :: cosE, senE
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: comp
    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: proj

    ! NELE: variável que armazena número de elementos
    ! NNOS: variável que armazena número de nós
    ! CELE: matriz que armazena índices dos nós iniciais e finais de cada elemento [nói, nóf]
    ! CNOS: posições x e y de cada nó [pos_x, pos_y]
    ! ROTS: matriz 3D que armazena a matriz de rotação local de cada elemento

    ! i: variável para ser iterada
    ! cosE: variável auxiliar para montagem da matriz de rotação local
    ! senE: variável auxiliar para montagem da matriz de rotação local
    ! proj: matriz que armazena projeções de cada elemento sobre o eixo x e y em m [proj_x, proj_y]
    ! comp: vetor que armazena comprimento de cada elemento em m

    ! Rotina para alocar vetores e matrizes
    ALLOCATE(comp(NELE))
    ALLOCATE(proj(NELE,2))

    ROTS = 0.0                      ! Iniciar matriz 3D com zeros

    DO i=1, NELE
        proj(i,1) = CNOS(CELE(i,2), 1) - CNOS(CELE(i,1), 1)         ! Projeção sobre o eixo x
        proj(i,2) = CNOS(CELE(i,2), 2) - CNOS(CELE(i,1), 2)         ! Projeção sobre o eixo y
        comp(i) = ((proj(i,1) ** 2) + (proj(i,2) ** 2)) ** 0.5      ! Comprimento

        cosE = proj(i,1) / comp(i)              ! Cosseno (projeção x sobre comprimento)
        senE = proj(i,2) / comp(i)              ! Seno (projeção y sobre comprimento)

        ! Início da matriz de rotação do elemento i

        ! Linha 1
        ROTS(1,1,i) = cosE
        ROTS(1,2,i) = senE

        ! Linha 2
        ROTS(2,1,i) = - senE
        ROTS(2,2,i) = cosE

        ! Linha 3
        ROTS(3,3,i) = 1

        ! Linha 4
        ROTS(4,4,i) = cosE
        ROTS(4,5,i) = senE

        ! Linha 5
        ROTS(5,4,i) = - senE
        ROTS(5,5,i) = cosE

        ! Linha 6
        ROTS(6,6,i) = 1

        ! Fim da matriz de rotação do elemento i
    END DO
END SUBROUTINE MMROT