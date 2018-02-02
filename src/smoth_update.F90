subroutine smoth_update()
  use openarray
  use variables
  use config
  implicit none
  
  q2b = q2+0.5d0*smoth*(q2f+q2b-2.d0*q2)
  q2 = q2f
  q2lb = q2l+0.5d0*smoth*(q2lf+q2lb-2.d0*q2l)
  q2l = q2lf
end subroutine smoth_update


subroutine smoth_update1(q2f_, u_, ub_)
  use openarray
  use variables
  use config
  implicit none
  type(array), intent(inout) :: q2f_, u_, ub_
  
  ub_ = u_+0.5d0*smoth*(q2f_+ub_-2.d0*u_)
  u_ = q2f_

end subroutine

subroutine fuck_you()
  use openarray
  use variables
  use config
  implicit none
  type(array) :: fk1, fk2, fk
  
    fk1 = ones(7,5,6,1,2)
    fk2 = ones(7,5,6,1,2)
    fk = fk1 + fk2

end subroutine

