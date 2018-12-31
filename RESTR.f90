SUBROUTINE RESTR(NNOS, REST, RIGG, MASG, AMOG, ESFG, CARM)
    ! Subrotina RESTR: Inserir RESTRições na matriz de rigidez global e vetor de esforços globais
    IMPLICIT NONE

    INTEGER, INTENT(in) :: NNOS
    LOGICAL, DIMENSION(3*NNOS), INTENT(in) :: REST
    REAL(kind=8), DIMENSION(3*NNOS) :: ESFG, CARM
    REAL(kind=8), DIMENSION(3*NNOS,3*NNOS) :: RIGG, MASG, AMOG

    INTEGER :: i, j

    ! NNOS: variável que armazena número de nós
    ! REST: vetor de valores logicos que significam restrição quando .TRUE.
    ! ESFG: vetor de esforços globais na estrutura
    ! RIGG: matriz de rigidez global da estrutura
    ! MASG: matriz de massa global da estrutura

    ! i: variável para iterar
    ! j: variável para iterar

    ! Inclusão das restrições
    DO i=1, 3*NNOS
        IF (REST(i)) THEN
            ESFG(i) = 0
            CARM(i) = 0
            DO j=1, 3*nnos
                IF (i == j) THEN
                    RIGG(i,i) = 1
                    MASG(i,i) = 1
                    AMOG(i,i) = 1
                ELSE
                    RIGG(i,j) = 0
                    RIGG(j,i) = 0
                    MASG(i,j) = 0
                    MASG(j,i) = 0
                    AMOG(i,j) = 0
                    AMOG(j,i) = 0
                ENDIF
            ENDDO
        ENDIF
    ENDDO

END SUBROUTINE RESTR