      Subroutine Grow1(Weight,U1,Xchain,Ychain,Zchain
     &    ,Array,Arr,Ar,Lovr,Lold,Rxx,Ryy,Rzz,P)



      Implicit None

      
      Integer, Parameter :: N = 10000,Nbead=6
      Logical Lovr,Lold
      Integer Ibead,Ran,Itrial,Ichoice,Itarget
      
       Double Precision Array(N),Arr(N),Ar(N)
     $  ,Xchain(Nbead),Ychain(Nbead),
     $   Zchain(Nbead),Weight,U1
     $  ,Dx,Dy,Dz,Ran_Uniform,R2,Nleft,pi,P(Nbead)
     $   ,Rxx(Nbead),Ryy(Nbead),Rzz(Nbead),Psum



       
       
       Xchain(1) = 0.0d0
       Ychain(1) = 0.0d0
       Zchain(1) = 0.0d0
       Xchain(Nbead) = 20.0d0
       Ychain(Nbead) = 0.0d0
       Zchain(Nbead) = 0.0d0
         Lovr = .False.
      Do Ibead = 2,Nbead-1
         Psum = 0.0d0
        Xchain(Ibead) = 0.0d0
        Ychain(Ibead) = 0.0d0
        Zchain(Ibead) = 0.0d0
      
          Dx = Xchain(Nbead) - Xchain(1)
          Dy = Ychain(Nbead) - Ychain(1)
          Dz = Zchain(Nbead) - Zchain(1)
           Dx = Dx/(Nbead-1)
           Dy = Dy/(Nbead-1)
           Dz = Dz/(Nbead-1)
          
          
          Xchain(Ibead) = Xchain(Ibead-1) + Dx
          Ychain(Ibead) = Ychain(Ibead-1) + Dy
          Zchain(Ibead) = Zchain(Ibead-1) + Dz
         Enddo

         
         End
