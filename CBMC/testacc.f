      Subroutine Testacc(Factor,Laccept)
      Implicit None

C     Test Acceptance

      Double Precision Ran_Uniform,Factor
      Logical Laccept
      
      If(Factor.Gt.1.0d0) Then

         Laccept = .True.

      Elseif(Factor.Lt.1.0d-9) Then

         Laccept = .False.

      Else

         If(Ran_Uniform().Lt.Factor) Then
            Laccept = .True.
         Else
            Laccept = .False.
         Endif

      Endif

      Return
      End
