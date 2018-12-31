cmake_minimum_required(VERSION 3.13)
project(trabalho05 Fortran)


enable_language(Fortran)
add_executable(
        trabalho05
        main.f90
        MMRIGL.f90
        MMROT.f90
        MMRIGG.f90
        MVESFL.f90
        MVESFG.f90
        RESTR.f90
        MVDESL.f90
        CVREAL.f90
        MMMASL.f90
        MMMASG.f90
        SPTRAN.f90
        CDELTAT.f90
        solver.f90)