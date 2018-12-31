SUBROUTINE MVESFL(NELE, NNOS, CELE, CNOS, CARE, ESFL)
    ! Subrotina MVESFL: Montar Vetor de ESForços Locais de cada elemento
    IMPLICIT NONE

    INTEGER, INTENT(in) :: NELE, NNOS
    INTEGER, DIMENSION(NELE,2), INTENT(in) :: CELE
    REAL(kind=8), DIMENSION(NNOS,2), INTENT(in) :: CNOS
    REAL(kind=8), DIMENSION(NELE,4), INTENT(in) :: CARE
    REAL(kind=8), DIMENSION(6,NELE), INTENT(out) :: ESFL

    INTEGER :: i
    REAL(kind=8), DIMENSION(:), ALLOCATABLE :: comp
    REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: proj

    ! NELE: variável que armazena número de elementos
    ! NNOS: variável que armazena número de nós
    ! CELE: matriz que armazena índices dos nós iniciais e finais de cada elemento [nói,nóf]
    ! CNOS: posições x e y de cada nó [pos_x,pos_y]
    ! CARE: matriz que armazena carregamentos locais de cada elemento em kN/m [Nli,Nlf,Qli,Qlf]
    ! ESFL: matriz 2D que armazena vetor de esforço local em cada elemento [6,num_elem]

    ! i: variável para ser iterada
    ! proj: matriz que armazena projeções de cada elemento sobre o eixo x e y em m [proj_x,proj_y]
    ! comp: vetor que armazena comprimento de cada elemento em m

    ! Rotina para alocar vetores e matrizes
    ALLOCATE(comp(NELE))
    ALLOCATE(proj(NELE,2))

    DO i=1, NELE
        proj(i,1) = CNOS(CELE(i,2), 1) - CNOS(CELE(i,1), 1)         ! Projeção sobre o eixo x
        proj(i,2) = CNOS(CELE(i,2), 2) - CNOS(CELE(i,1), 2)         ! Projeção sobre o eixo y
        comp(i) = ((proj(i,1) ** 2) + (proj(i,2) ** 2)) ** 0.5      ! Comprimento

        ! Montar matriz de esforços locais para cada elemento
        ESFL(1,i) = (CARE(i,1) / 3 + CARE(i,2) / 6) * comp(i)               ! (Ni/3 + Nf/6)*L
        ESFL(2,i) = (7 * CARE(i,3) / 20 + 3 * CARE(i,4) / 20) * comp(i)     ! (7*Qi/20 + 3Qf/20)*L
        ESFL(3,i) = (CARE(i,3) / 20 + CARE(i,4) / 30) * comp(i) ** 2        ! (Qi/20 + Qf/30)*L^2
        ESFL(4,i) = (CARE(i,1) / 6 + CARE(i,2) / 3) * comp(i)               ! (Ni/6 + Nf/3)*L
        ESFL(5,i) = (3 * CARE(i,3) / 20 + 7 * CARE(i,4) / 20) * comp(i)     ! (3*Qi/20 + 7Qf/20)*L
        ESFL(6,i) = (- CARE(i,3) / 30 - CARE(i,4) / 20) * comp(i) ** 2      ! (-Qi/30 - Qf/20)*L^2
    END DO
END SUBROUTINE MVESFL