       Program Mainprogram
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Main program                                                    C
C     CBMC simulation of single chain nanoparticles                   C
C     Bernardo Oyarzun - 2017                                         C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      Implicit None
     
      
      Integer, Parameter :: N = 10000, Nbead = 11,
     $      Ncycle = 10000, Nrun = 1
      Integer I,Ibead,Accept,Icycle,Irun
      Logical Lovr,Lold,Laccept
      Double Precision Array(N),Arr(N),Ar(N)
     $  ,Rxx(Nbead),Ryy(Nbead),Rzz(Nbead),Xchain(Nbead)
     $  ,Weight,U1,Ychain(Nbead),Zchain(Nbead)
     $  ,Dx,Dy,Dz,R2,U,P(Nbead),Uold,Rxx_New(Nbead)
     $  ,Ryy_New(Nbead),Rzz_New(Nbead),Wnew,Wold
     $  ,Pold,Unew,Pnew,U_new,U_old,R2new,R2old
    
    
     

      Call Gaussian(Array,Arr,Ar)
1     Continue
      
      Call Grow1(Weight,U1,Xchain,Ychain,Zchain,
     &   Array,Arr,Ar,Lovr,Lold,Rxx,Ryy,Rzz,P)
       
      If(Lovr) Goto 1
      Do Ibead=1,Nbead
      Rxx(Ibead) = Xchain(Ibead)
      Ryy(Ibead) = Ychain(Ibead)
      Rzz(Ibead) = Zchain(Ibead)
      Enddo
      Open(24,File='Xcord2.txt',Status='Unknown')
      Open(25,File='Xcord5.txt',Status='Unknown')
      Open(26,File='Xcord9.txt',Status='Unknown')
      Open(22,File='New_coord.txt',Status="Unknown")

      Accept = 0
C    Growing a new configuration
      Do Irun = 1, Nrun
       Accept = 0
      Do Icycle=1,Ncycle
      Lold = .False.
     
       Pnew = 1.0d0

      Call Grow(Weight,U1,Xchain,Ychain,Zchain,
     &   Array,Arr,Ar,Lovr,Lold,Rxx,Ryy,Rzz,Pnew)

      
       Do Ibead=1,Nbead
        Rxx_New(Ibead) = Xchain(Ibead)
        Ryy_New(Ibead) = Ychain(Ibead)
        Rzz_New(Ibead) = Zchain(Ibead)
        Enddo
           
         
        
        Dx = Rxx_New(Nbead) - Rxx_New(Nbead-1)
        Dy = Ryy_New(Nbead) - Ryy_New(Nbead-1)
        Dz = Rzz_New(Nbead) - Rzz_New(Nbead-1)
        R2 = Dx**2 + Dy**2 + Dz**2
         R2new = R2
        Call energy(U,Lovr,R2)
         U_new = U
        
         
          
        
       Wnew = Weight

C  growing an old configuration

      Lold = .True.
      Pold = 1.0d0
      
      Call Grow(Weight,U1,Xchain,Ychain,Zchain,
     &   Array,Arr,Ar,Lovr,Lold,Rxx,Ryy,Rzz,Pold)

      
      Do Ibead=1,Nbead
      Rxx(Ibead) = Xchain(Ibead)
      Ryy(Ibead) = Ychain(Ibead)
      Rzz(Ibead) = Zchain(Ibead)
      Enddo
      Uold = U1
      Wold = Weight
      
     
      Dx = Rxx(Nbead) - Rxx(Nbead-1)
      Dy = Ryy(Nbead) - Ryy(Nbead-1)
      Dz = Rzz(Nbead) - Rzz(Nbead-1)
      R2 = Dx**2 + Dy**2 + Dz**2
       R2old = R2
      Call energy(U,Lovr,R2)
       U_old = U
     
     
      Call Testacc((Wnew*exp(-U_new+U_old)*Pold)
     &          /(Wold*Pnew),Laccept)
      
      If(Laccept) Then
        Accept = Accept + 1
       
       Do I = 1,Nbead
        Rxx(I) = Rxx_New(I)
        Ryy(I) = Ryy_New(I)
        Rzz(I) = Rzz_New(I)
       Enddo
       
        
       Endif
       
         
          Write(26,*) Rxx(5)
          Write(24,*) Rxx(2)
         
          
          Write(25,*) Rxx(3)
         
      
        
      
       
       Enddo
        Open(23,File='Accept.txt',Status="Unknown")
       Write(23,*) Accept
        Enddo
        Close(24)
        Close(23)
        Close(24)
       End



     
