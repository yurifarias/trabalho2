PROGRAM trabalho
    IMPLICIT NONE
    
    CHARACTER*20 titulo
    CHARACTER*96 detalhes
    INTEGER arq1, arq2, numnos, numele, i, j, k, l
    INTEGER noi[ALLOCATABLE](:) , nof[ALLOCATABLE](:)
    REAL*8 val1, val2, val3, val4, val5, cosE, senE
    REAL*8 mate[ALLOCATABLE](:), larg[ALLOCATABLE](:), altu[ALLOCATABLE](:), area[ALLOCATABLE](:), iner[ALLOCATABLE](:)
    REAL*8 posx[ALLOCATABLE](:), posy[ALLOCATABLE](:), compx[ALLOCATABLE](:), compy[ALLOCATABLE](:), comp[ALLOCATABLE](:)
    DOUBLE PRECISION mrloc(6,6), mrot(6,6), mrott(6,6), mrelglo(6,6), maux(6,6)
    
    arq1 = 1
    arq2 = 2
    
    ! Abrir arquivos de entrada e saida de dados
    OPEN(arq1, FILE='entrada.txt', STATUS='unknown')
    OPEN(arq2, FILE='saida.txt', STATUS='unknown')
    
    ! Ler titulo do arquivo de entrada
    READ(arq1,'(A)') titulo
    WRITE(arq2,*) titulo
    WRITE(arq2,*)
    
    READ(arq1,*)
    READ(arq1,*)
    
    ! Ler numero de elementos e definir tamanho de vetores
    READ(arq1,*) numele   
    READ(arq1,*)
    
    READ(arq1,'(A)') detalhes
    WRITE(arq2,'(A)') detalhes
    
    ALLOCATE(mate(numele), larg(numele), altu(numele), noi(numele), nof(numele))
    
    ! Ler dados de vetores
    DO i=1,numele
        READ(arq1,*) j, noi(j), nof(j), mate(j), larg(j), altu(j)
        WRITE(arq2,12) j, noi(j), nof(j), mate(j), larg(j), altu(j)
    ENDDO
    
    READ(arq1,*)
    READ(arq1,*)
    WRITE(arq2,*)
    
    ! Ler numero de nos e definir tamanho de vetores
    READ(arq1,*) numnos
    
    READ(arq1,*)
    READ(arq1,'(A)') detalhes
    WRITE(arq2,'(A)') detalhes
    
    ALLOCATE(posx(numnos), posy(numnos))
    
    ! Ler dados de vetores
    DO i=1, numnos
        READ(arq1,*) j, posx(j), posy(j) 
        WRITE(arq2,11) j, posx(j), posy(j)
    ENDDO
    
    ALLOCATE(area(numele), iner(numele), compx(numele), compy(numele), comp(numele))
    
    ! Calcular area, inercia e comprimento pra cada elemento, assim como sua matriz de rigidez e de rotacao
    DO i=1, numele
        area(i) = altu(i) * larg(i)
        iner(i) = altu(i) * (larg(i) ** 3) / 12
        compx(i) = posx(nof(i)) - posx(noi(i))
        compy(i) = posy(nof(i)) - posy(noi(i))
        comp(i) = ((compx(i) ** 2) + (compy(i) ** 2)) ** 0.5
        
        val1 = mate(i) * area(i) / comp(i)
        val2 = 12 * mate(i) * iner(i) / (comp(i) ** 3)
        val3 = 6 * mate(i) * iner(i) / (comp(i) ** 2)
        val4 = 4 * mate(i) * iner(i) / comp(i)
        val5 = 2 * mate(i) * iner(i) / comp(i)
        
        cosE = (posx(nof(i)) - posx(noi(i))) / comp(i)
        senE = (posy(nof(i)) - posy(noi(i))) / comp(i)
        
        mrloc(1,1) = val1
        mrloc(1,2) = 0
        mrloc(1,3) = 0
        mrloc(1,4) = - val1
        mrloc(1,5) = 0
        mrloc(1,6) = 0
        
        mrloc(2,1) = 0
        mrloc(2,2) = val2
        mrloc(2,3) = val3
        mrloc(2,4) = 0
        mrloc(2,5) = - val2
        mrloc(2,6) = val3
        
        mrloc(3,1) = 0
        mrloc(3,2) = val3
        mrloc(3,3) = val4
        mrloc(3,4) = 0
        mrloc(3,5) = - val3
        mrloc(3,6) = val5
        
        mrloc(4,1) = - val1
        mrloc(4,2) = 0
        mrloc(4,3) = 0
        mrloc(4,4) = val1
        mrloc(4,5) = 0
        mrloc(4,6) = 0
        
        mrloc(5,1) = 0
        mrloc(5,2) = - val2
        mrloc(5,3) = - val3
        mrloc(5,4) = 0
        mrloc(5,5) = val2
        mrloc(5,6) = - val3
        
        mrloc(6,1) = 0
        mrloc(6,2) = val3
        mrloc(6,3) = val5
        mrloc(6,4) = 0
        mrloc(6,5) = - val3
        mrloc(6,6) = val4
        
        mrot(1,1) = cosE
        mrot(1,2) = senE
        mrot(1,3) = 0
        mrot(1,4) = 0
        mrot(1,5) = 0
        mrot(1,6) = 0
        
        mrot(2,1) = - senE
        mrot(2,2) = cosE
        mrot(2,3) = 0
        mrot(2,4) = 0
        mrot(2,5) = 0
        mrot(2,6) = 0
        
        mrot(3,1) = 0
        mrot(3,2) = 0
        mrot(3,3) = 1
        mrot(3,4) = 0
        mrot(3,5) = 0
        mrot(3,6) = 0
        
        mrot(4,1) = 0
        mrot(4,2) = 0
        mrot(4,3) = 0
        mrot(4,4) = cosE
        mrot(4,5) = senE
        mrot(4,6) = 0
        
        mrot(5,1) = 0
        mrot(5,2) = 0
        mrot(5,3) = 0
        mrot(5,4) = - senE
        mrot(5,5) = cosE
        mrot(5,6) = 0
        
        mrot(6,1) = 0
        mrot(6,2) = 0
        mrot(6,3) = 0
        mrot(6,4) = 0
        mrot(6,5) = 0
        mrot(6,6) = 1
        
        DO k=1, 6
            DO j=1, 6
                mrott(k,j) = mrot(j,k)
            ENDDO
        ENDDO
        
        WRITE(arq2,*)
        WRITE(arq2,*)
        
        WRITE(arq2,18) 'ELEMENTO', i
        WRITE(arq2,*)
        
        WRITE(arq2,'(A)') 'COMPRIMENTO (m) COSSENO         SENO            AREA (m2)       MINERCIA (m4)   MELASTICIDADE (N/m2)'
        WRITE(arq2,*)
        WRITE(arq2,17) comp(i), cosE, senE, area(i), iner(i), mate(i)
        WRITE(arq2,*)
        
        WRITE(arq2,14) 'MATRIZ DE RIGIDEZ LOCAL DO ELEMENTO ', i
        WRITE(arq2,*)
        DO k=1, 6
            WRITE(arq2, 13) mrloc(k,1), mrloc(k,2), mrloc(k,3), mrloc(k,4), mrloc(k,5), mrloc(k,6)
        ENDDO
        WRITE(arq2,*)
        
        WRITE(arq2,15) 'MATRIZ DE ROTACAO DO ELEMENTO ', i
        WRITE(arq2,*)
        DO k=1, 6
            WRITE(arq2, 13) mrot(k,1), mrot(k,2), mrot(k,3), mrot(k,4), mrot(k,5), mrot(k,6)
        ENDDO
        WRITE(arq2,*)
        
        WRITE(arq2,16) 'MATRIZ DE ROTACAO TRANSPOSTA DO ELEMENTO ', i
        WRITE(arq2,*)
        DO k=1, 6
            WRITE(arq2, 13) mrott(k,1), mrott(k,2), mrott(k,3), mrott(k,4), mrott(k,5), mrott(k,6)
        ENDDO
        WRITE(arq2,*)
        
        DO j=1, 6
            DO k=1, 6
                maux(j,k) = 0
                DO l=1,6
                    maux(j,k) = maux(j,k) + mrott(j,l) * mrloc(l,k)
                ENDDO
            ENDDO
        ENDDO
        
        DO j=1, 6
            DO k=1, 6
                mrelglo(j,k) = 0
                DO l=1,6
                    mrelglo(j,k) = mrelglo(j,k) + maux(j,l) * mrot(l,k)
                ENDDO
            ENDDO
        ENDDO
        
        WRITE(arq2,14) 'MATRIZ DE RIGIDEZ GLOBAL DO ELEMENTO ', i
        WRITE(arq2,*)
        DO k=1, 6
            WRITE(arq2, 13) mrelglo(k,1), mrelglo(k,2), mrelglo(k,3), mrelglo(k,4), mrelglo(k,5), mrelglo(k,6)
        ENDDO
        
    ENDDO
    
    CLOSE(arq1)
    CLOSE(arq2)
    
10 FORMAT(3x,3a8)    
11 FORMAT(1i16,2f16.3)
12 FORMAT(3i16,3f16.3)
13 FORMAT(6f16.3)
14 FORMAT(a37,i2)
15 FORMAT(a30,i2) 
16 FORMAT(a41,i2)
17 FORMAT(es16.3,2f16.3,3es16.3)
18 FORMAT(a8, i2)   
    
END PROGRAM trabalho