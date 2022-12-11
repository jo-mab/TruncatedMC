       Program Mainprogram
      


      Implicit None
     
      
      Integer, Parameter :: N = 10000, Nbead = 6,
     $      Ncycle = 10000, Nrun = 1
      Integer I,Ibead,Accept,Icycle,Irun
      Logical Lovr,Lold,Laccept,Accept1(Ncycle)
      Double Precision Array(N),Arr(N),Ar(N)
     $  ,Rxx(Nbead),Ryy(Nbead),Rzz(Nbead),Xchain(Nbead)
     $  ,Weight,U1,Ychain(Nbead),Zchain(Nbead)
     $  ,Dx,Dy,Dz,R2,U,P(Nbead),Uold,Rxx_New(Nbead)
     $  ,Ryy_New(Nbead),Rzz_New(Nbead),Wnew,Wold
     $  ,Pold,Unew,Pnew,U_new,U_old,R2new,R2old,
     $   Waverage,Wnew1(Ncycle),standard_error,sigma,
     $   Wold1(Ncycle), Waverage_Equilibrium,
     $   Waverage_Nequilibrium,D,Wnew2(Ncycle),
     $   Wcheck,Pnew1(Ncycle), Pold1(Ncycle)
     
    
    
     

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
      Open(24,File='Accepted.txt',Status='Unknown')
      Open(25,File='Waverage.txt',Status='Unknown')
      Open(26,File='W_error.txt',Status='Unknown')
      
      Open(27,File='Wcycle_Equilibrium.txt',Status="Unknown")
      Open(28,File='Wcycle_NEquilibrium.txt',Status="Unknown")
      Dx=Rxx(Nbead)-Rxx(1)
      Dy=Ryy(Nbead)-Ryy(1)
      Dz=Rzz(Nbead)-Rzz(1)
      D = Dx*Dx + Dy*Dy +Dz*Dz
      Write(*,*) D
      Accept = 0
      Waverage=0
C    Growing a new configuration
      Do Irun = 1, Nrun
       Accept = 0
      Do Icycle=1,Ncycle
      Lold = .False.
      Accept1(Icycle) = .False.
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
      write(25,*) Icycle, Wnew, Wold
      Wnew1(Icycle)=Wnew*exp(-U_new)
      Wold1(Icycle) = Wold*exp(-U_old)
      
      Pnew1(Icycle) = Pnew
      Pold1(Icycle) = Pold
     
      Call Testacc((Wnew*exp(-U_new+U_old)*Pold)
     &          /(Wold*Pnew),Laccept)
      
      If(Laccept) Then
        Accept = Accept + 1
        Accept1(Icycle) = .True.
        write(24,*) Icycle, Wnew, Wold
        
       Do I = 1,Nbead
        Rxx(I) = Rxx_New(I)
        Ryy(I) = Ryy_New(I)
        Rzz(I) = Rzz_New(I)
       Enddo
       
        
       Endif
       
         
         
         
      
        
      
       
       Enddo
        Open(23,File='Accept.txt',Status="Unknown")
       Write(23,*) Accept
        Enddo
        Waverage_Equilibrium=0
        Waverage_Nequilibrium=0
        Do Icycle=1,Ncycle
           Waverage_Nequilibrium = Waverage_Nequilibrium+
     &             Wnew1(Icycle)
           Wcheck = Wnew1(Icycle)/Pnew1(Icycle)
           Write(28,*) Wcheck
           If(Accept1(Icycle)) Then
               Waverage_Equilibrium=Waverage_Equilibrium
     &              +Wnew1(Icycle)
               Wcheck = Wnew2(Icycle)/Pnew1(Icycle)
               Write(27,*) Wcheck
           Else
               Waverage_Equilibrium=Waverage_Equilibrium
     &               +Wold1(Icycle)
               Wcheck = Wold1(Icycle)/Pold1(Icycle)
               Write(27,*) Wcheck
           Endif

           
        Enddo
        Waverage_Equilibrium = Waverage_Equilibrium/Ncycle
        Waverage_NEquilibrium = Waverage_NEquilibrium/Ncycle
        
        
       
        Write(26,*) -log(Waverage_Equilibrium),
     &              -log(Waverage_Nequilibrium)
        
        Close(25)
        Close(23)
        Close(24)
        Close(26)
       End



     
