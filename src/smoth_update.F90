subroutine smoth_update(uf_, vf_, u_, v_, ub_, vb_)
  use openarray
  use variables
  use config
  implicit none
  type(array), intent(inout) :: uf_, vf_, u_, v_, ub_, vb_
  
  ub_%data = u_%data+0.5d0*smoth*(uf_%data+ub_%data-2.d0*u_%data)
  
  u_%data = uf_%data
  
  vb_%data = v_%data+0.5d0*smoth*(vf_%data+vb_%data-2.d0*v_%data)
  
  v_%data = vf_%data

end subroutine smoth_update


subroutine smoth_update1(uf_, u_, ub_)
  use openarray
  use variables
  use config
  implicit none
  type(array), intent(inout) :: uf_, u_, ub_
  
  ub_%data = u_%data+0.5d0*smoth*(uf_%data+ub_%data-2.d0*u_%data)
  u_%data = uf_%data

end subroutine
