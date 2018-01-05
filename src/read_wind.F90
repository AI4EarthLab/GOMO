
subroutine read_var(var, var_name, prefix,  iint, iw)
  use variables
  use config
  use openarray
  implicit none
  integer :: iint
  integer :: iw
  type(array), intent(inout) :: var
  character(len=1000) :: fnc
  character(len=*) :: file_prefix, var_name
  
  fnc = trim(in_path)//trim(i2s(problem)) &
       //"."//prefix//"."&
       //trim(i2s(int(iint/iw), '(I0.4)'))//".nc"

  var = load(fnc, var_name)
  
end subroutine
