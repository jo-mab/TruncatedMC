      Subroutine Gaussian (Array,Arr,Ar)
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Main program                                                    C
C     CBMC simulation of single chain nanoparticles                   C
C     Bernardo Oyarzun - 2017                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Implicit None
     
      
      
       Integer, Parameter :: N = 10000
       Integer i
       Double Precision Array(N),Pi,Temp, Mean,SD,Arr(N),Ar(N)

      Pi = 4.0*Atan(1.0)
      Mean = 0
      SD = 1.0
      Call Random_Number(Array)
       
      Do i=1,N-1,2
      Temp = SD * SQRT(-2.0*LOG(Array(i))) * COS(2*Pi*Array(i+1)) + Mean
       Array(i+1) = SD * SQRT(-2.0*LOG(Array(i))) * SIN(2*pi*Array(i+1))
     &     + Mean
       Array(i) = Temp

      Enddo
      
      Call Random_Number(Arr)
      Do i=1,N-1,2
      Temp = SD * SQRT(-2.0*LOG(Arr(i))) * COS(2*Pi*Arr(i+1))+Mean
      Arr(i+1) = SD * SQRT(-2.0*LOG(Arr(i))) * SIN(2*pi*Arr(i+1))
     &     + Mean
      
      Arr(i) = Temp
      Enddo
 
      Call Random_Number(Ar)
       Do i=1,N-1,2
        Temp = SD * SQRT(-2.0*LOG(Ar(i))) * COS(2*Pi*Ar(i+1))
     &    + Mean
        Ar(i+1) = SD * SQRT(-2.0*LOG(Ar(i))) * SIN(2*pi*Ar(i+1))
     &      + Mean
        Ar(i) = Temp
       Enddo
      
            
       
       End




