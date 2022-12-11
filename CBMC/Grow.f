      Subroutine Grow(Weight,U1,Xchain,Ychain,Zchain
     &    ,Array,Arr,Ar,Lovr,Lold,Rxx,Ryy,Rzz,Pguiding)



      Implicit None

      
      Integer, Parameter :: N = 10000,Nbead=6, Ktrial = 1000
      Logical Lovr,Loverlap(Ktrial),Lold
      Integer Ibead,Ran,Itrial,Ichoice,Itarget,I,J,Ichoice1
      
       Double Precision Array(N),Arr(N),Ar(N)
     $  ,Xchain(Nbead),Ychain(Nbead),
     $   Zchain(Nbead),Weight,U1
     $  ,Xtrial(Ktrial),Ytrial(Ktrial),
     $   Ztrial(Ktrial),W(Ktrial),U,Psum,Pmax
     $   ,Utrial(Ktrial),Sumw,Ptrial(Ktrial)
     $   ,Beta,Dx,Dy,Dz,Ran_Uniform,R2,Nleft,pi,P(Nbead)
     $   ,Rxx(Nbead),Ryy(Nbead),Rzz(Nbead),Rsquare,Eval
     $   ,Dx1,Dy1,Dz1,Wsum,Psum1,y,Pguide(Ktrial),Nright,Rsqu,
     $    Pguiding,Number,Sum,s,s1,Weight1,Wsum1


       
       pi = 4.0*Atan(1.0)
       Lovr = .False.
       Weight = 1.0d0
       Weight1=1.0d0
       U1 = 0.0d0
       Pguiding = 1.0d0
       Beta = 1.0d0
       Xchain(1) = 0.0d0
       Ychain(1) = 0.0d0
       Zchain(1) = 0.0d0
       Xchain(Nbead) = 20.0d0
       Ychain(Nbead) = 0.0d0
       Zchain(Nbead) = 0.0d0
       
      Do Ibead = 2,Nbead-1
         Psum = 0.0d0
        Xchain(Ibead) = 0.0d0
        Ychain(Ibead) = 0.0d0
        Zchain(Ibead) = 0.0d0
    
        Do Itrial = 1, Ktrial
         Utrial(Itrial) = 0.0d0
         U = 0.0d0
         If(Lold.And.(Itrial.eq.1)) Then
           Xtrial(Itrial) = Rxx(Ibead)
           Ytrial(Itrial) = Ryy(Ibead)
           Ztrial(Itrial) = Rzz(Ibead)
           
          Else

         
          call random_number(s)
          call random_number(s1)
          Dx = sqrt(-2.0*log(1-s)) * Cos(2*Pi*s1)
          call random_number(s)
          call random_number(s1)
           Dy = sqrt(-2.0*log(1-s)) * Cos(2*Pi*s1)
          call random_number(s)
          call random_number(s1)
           Dz = sqrt(-2.0*log(1-s)) * Cos(2*Pi*s1)
          
          Xtrial(Itrial) = Xchain(Ibead-1) + Dx
          Ytrial(Itrial) = Ychain(Ibead-1) + Dy
          Ztrial(Itrial) = Zchain(Ibead-1) + Dz
         
         Endif
         
         
            Dx1 = Xtrial(Itrial) - Xchain(Nbead)
            Dy1 = Ytrial(Itrial) - Ychain(Nbead)
            Dz1 = Ztrial(Itrial) - Zchain(Nbead)
          
            Rsquare = Dx1**2 + Dy1**2 + Dz1**2
              Nleft = Nbead - Ibead
             
            
             Number = (2*pi*Nleft)**(1.5)
             
            Ptrial(Itrial) = exp(-1*Rsquare/(2*Nleft))
            Pguide(Itrial) = Ptrial(Itrial)
            
           Dx1 = Xtrial(Itrial) - Xchain(Ibead-1)
           Dy1 = Ytrial(Itrial) - Ychain(Ibead-1)
           Dz1 = Ztrial(Itrial) - Zchain(Ibead-1)
          
           R2 = Dx1*Dx1 + Dy1*Dy1 + Dz1*Dz1
           Call energy(U,Lovr,R2)
           Utrial(Itrial) = U
           
            
         Enddo
           
         
          
         
C     Select trial configuration
         Wsum = 0.0d0
         Wsum1=0.0d0
        Do Itrial =1,Ktrial
        
        
         Wsum = Wsum + Pguide(Itrial)
         
        Enddo
       
         
          Sum = 0.0d0
         Do Itrial = 1, Ktrial
          Ptrial(Itrial) = Pguide(Itrial)/Wsum
           Sum = Sum+Ptrial(Itrial)
          
         Enddo
         
        
         y = Ran_Uniform()
         
         
        
         
          Ichoice = 0
         Do Itrial=1,Ktrial
           Psum = 0.0d0
           Psum1 = 0.0d0
          Do I=1,Itrial-1
           If(I.eq.0) Then
            Psum = 0.0d0
          Else
           Psum = Psum + Ptrial(I)
           
          Endif
          Enddo
          Do I=1,Itrial
           Psum1= Psum1+Ptrial(I)
          
          Enddo
          
          If((y.Gt.Psum).And.(y.Lt.Psum1)) Then
            Ichoice = Itrial
            EXIT
          Endif
        Enddo
         If(Lold) Ichoice = 1
        
        
       Wsum = Wsum/ktrial
        Weight = Weight*Wsum
       
        
        
        
        Pguiding = Pguiding*Pguide(Ichoice)
        
       
        Xchain(Ibead) = Xtrial(Ichoice)
        Ychain(Ibead) = Ytrial(Ichoice)
        Zchain(Ibead) = Ztrial(Ichoice)
         
       
        
        Enddo
         End
