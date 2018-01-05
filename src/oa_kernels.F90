#include "config.h"
module oa_kernels
  use oa_dict
  use oa_petsc
  use oa_type
  use oa_array
  use oa_buffer
  use oa_transfer
  use oa_expr

contains

  subroutine init_kernels(ierr)
    implicit none
    integer,intent(out) :: ierr

    !call dict_add(dict_kernels, 5319368613274967336_8, loc(kernel_5319368613274967336))    
    call dict_add(dict_kernels, 5319368613274967336_8, kernel_5319368613274967336)
    !call dict_add(dict_kernels, 4639977022766068446_8, loc(kernel_4639977022766068446))    
    call dict_add(dict_kernels, 4639977022766068446_8, kernel_4639977022766068446)
    !call dict_add(dict_kernels, 5104477621180273533_8, loc(kernel_5104477621180273533))    
    call dict_add(dict_kernels, 5104477621180273533_8, kernel_5104477621180273533)
    !call dict_add(dict_kernels, 3052836676995476825_8, loc(kernel_3052836676995476825))    
    call dict_add(dict_kernels, 3052836676995476825_8, kernel_3052836676995476825)
    !call dict_add(dict_kernels, 3254872525222050683_8, loc(kernel_3254872525222050683))    
    call dict_add(dict_kernels, 3254872525222050683_8, kernel_3254872525222050683)
    !call dict_add(dict_kernels, 4066808250212656659_8, loc(kernel_4066808250212656659))    
    call dict_add(dict_kernels, 4066808250212656659_8, kernel_4066808250212656659)
    !call dict_add(dict_kernels, 3262871250189963503_8, loc(kernel_3262871250189963503))    
    call dict_add(dict_kernels, 3262871250189963503_8, kernel_3262871250189963503)
    !call dict_add(dict_kernels, 3384383430676826436_8, loc(kernel_3384383430676826436))    
    call dict_add(dict_kernels, 3384383430676826436_8, kernel_3384383430676826436)
    !call dict_add(dict_kernels, 1919310359753992010_8, loc(kernel_1919310359753992010))    
    call dict_add(dict_kernels, 1919310359753992010_8, kernel_1919310359753992010)
    !call dict_add(dict_kernels, 5522439168916705770_8, loc(kernel_5522439168916705770))    
    call dict_add(dict_kernels, 5522439168916705770_8, kernel_5522439168916705770)
    !call dict_add(dict_kernels, 3692924362984551671_8, loc(kernel_3692924362984551671))    
    call dict_add(dict_kernels, 3692924362984551671_8, kernel_3692924362984551671)
    !call dict_add(dict_kernels, 7913083383874115848_8, loc(kernel_7913083383874115848))    
    call dict_add(dict_kernels, 7913083383874115848_8, kernel_7913083383874115848)
    !call dict_add(dict_kernels, 646399903993819770_8, loc(kernel_646399903993819770))    
    call dict_add(dict_kernels, 646399903993819770_8, kernel_646399903993819770)
    !call dict_add(dict_kernels, 7761457064372070790_8, loc(kernel_7761457064372070790))    
    call dict_add(dict_kernels, 7761457064372070790_8, kernel_7761457064372070790)
    !call dict_add(dict_kernels, 644868325008555321_8, loc(kernel_644868325008555321))    
    call dict_add(dict_kernels, 644868325008555321_8, kernel_644868325008555321)
    !call dict_add(dict_kernels, 646399903991447928_8, loc(kernel_646399903991447928))    
    call dict_add(dict_kernels, 646399903991447928_8, kernel_646399903991447928)
    !call dict_add(dict_kernels, 6766139683933094313_8, loc(kernel_6766139683933094313))    
    call dict_add(dict_kernels, 6766139683933094313_8, kernel_6766139683933094313)
    !call dict_add(dict_kernels, 108513903478386821_8, loc(kernel_108513903478386821))    
    call dict_add(dict_kernels, 108513903478386821_8, kernel_108513903478386821)
    !call dict_add(dict_kernels, 2381965364436835615_8, loc(kernel_2381965364436835615))    
    call dict_add(dict_kernels, 2381965364436835615_8, kernel_2381965364436835615)
    !call dict_add(dict_kernels, 5084845002459917026_8, loc(kernel_5084845002459917026))    
    call dict_add(dict_kernels, 5084845002459917026_8, kernel_5084845002459917026)
    !call dict_add(dict_kernels, 5404229555782944892_8, loc(kernel_5404229555782944892))    
    call dict_add(dict_kernels, 5404229555782944892_8, kernel_5404229555782944892)
    !call dict_add(dict_kernels, 644868325004997558_8, loc(kernel_644868325004997558))    
    call dict_add(dict_kernels, 644868325004997558_8, kernel_644868325004997558)
    !call dict_add(dict_kernels, 7381972720285064861_8, loc(kernel_7381972720285064861))    
    call dict_add(dict_kernels, 7381972720285064861_8, kernel_7381972720285064861)
    !call dict_add(dict_kernels, 8205386001274835749_8, loc(kernel_8205386001274835749))    
    call dict_add(dict_kernels, 8205386001274835749_8, kernel_8205386001274835749)
    !call dict_add(dict_kernels, 5507350668782427326_8, loc(kernel_5507350668782427326))    
    call dict_add(dict_kernels, 5507350668782427326_8, kernel_5507350668782427326)
    !call dict_add(dict_kernels, 3303401064826175212_8, loc(kernel_3303401064826175212))    
    call dict_add(dict_kernels, 3303401064826175212_8, kernel_3303401064826175212)
    !call dict_add(dict_kernels, 2281688410210494257_8, loc(kernel_2281688410210494257))    
    call dict_add(dict_kernels, 2281688410210494257_8, kernel_2281688410210494257)
    !call dict_add(dict_kernels, 259121922155312303_8, loc(kernel_259121922155312303))    
    call dict_add(dict_kernels, 259121922155312303_8, kernel_259121922155312303)
    !call dict_add(dict_kernels, 7595890524641150240_8, loc(kernel_7595890524641150240))    
    call dict_add(dict_kernels, 7595890524641150240_8, kernel_7595890524641150240)
    !call dict_add(dict_kernels, 6383578275131078730_8, loc(kernel_6383578275131078730))    
    call dict_add(dict_kernels, 6383578275131078730_8, kernel_6383578275131078730)
    !call dict_add(dict_kernels, 8877419281824194841_8, loc(kernel_8877419281824194841))    
    call dict_add(dict_kernels, 8877419281824194841_8, kernel_8877419281824194841)
    !call dict_add(dict_kernels, 1039351428668170720_8, loc(kernel_1039351428668170720))    
    call dict_add(dict_kernels, 1039351428668170720_8, kernel_1039351428668170720)
    !call dict_add(dict_kernels, 1279377163763675870_8, loc(kernel_1279377163763675870))    
    call dict_add(dict_kernels, 1279377163763675870_8, kernel_1279377163763675870)
    !call dict_add(dict_kernels, 8707219604264828824_8, loc(kernel_8707219604264828824))    
    call dict_add(dict_kernels, 8707219604264828824_8, kernel_8707219604264828824)
    !call dict_add(dict_kernels, 4905021281077427515_8, loc(kernel_4905021281077427515))    
    call dict_add(dict_kernels, 4905021281077427515_8, kernel_4905021281077427515)
    !call dict_add(dict_kernels, 1819205572515252931_8, loc(kernel_1819205572515252931))    
    call dict_add(dict_kernels, 1819205572515252931_8, kernel_1819205572515252931)
    !call dict_add(dict_kernels, 5940693024572846470_8, loc(kernel_5940693024572846470))    
    call dict_add(dict_kernels, 5940693024572846470_8, kernel_5940693024572846470)
    !call dict_add(dict_kernels, 3864818436564201008_8, loc(kernel_3864818436564201008))    
    call dict_add(dict_kernels, 3864818436564201008_8, kernel_3864818436564201008)
    !call dict_add(dict_kernels, 1820361344622479259_8, loc(kernel_1820361344622479259))    
    call dict_add(dict_kernels, 1820361344622479259_8, kernel_1820361344622479259)
    !call dict_add(dict_kernels, 2483251290183065470_8, loc(kernel_2483251290183065470))    
    call dict_add(dict_kernels, 2483251290183065470_8, kernel_2483251290183065470)
    !call dict_add(dict_kernels, 8238309252147448007_8, loc(kernel_8238309252147448007))    
    call dict_add(dict_kernels, 8238309252147448007_8, kernel_8238309252147448007)
    !call dict_add(dict_kernels, 445431549541144228_8, loc(kernel_445431549541144228))    
    call dict_add(dict_kernels, 445431549541144228_8, kernel_445431549541144228)
    !call dict_add(dict_kernels, 6348923055873632844_8, loc(kernel_6348923055873632844))    
    call dict_add(dict_kernels, 6348923055873632844_8, kernel_6348923055873632844)
    !call dict_add(dict_kernels, 4164740083217831399_8, loc(kernel_4164740083217831399))    
    call dict_add(dict_kernels, 4164740083217831399_8, kernel_4164740083217831399)
    !call dict_add(dict_kernels, 7821363202111267096_8, loc(kernel_7821363202111267096))    
    call dict_add(dict_kernels, 7821363202111267096_8, kernel_7821363202111267096)
    !call dict_add(dict_kernels, 7381297727881908921_8, loc(kernel_7381297727881908921))    
    call dict_add(dict_kernels, 7381297727881908921_8, kernel_7381297727881908921)
    !call dict_add(dict_kernels, 563110658223512454_8, loc(kernel_563110658223512454))    
    call dict_add(dict_kernels, 563110658223512454_8, kernel_563110658223512454)
    !call dict_add(dict_kernels, 3845637904882483465_8, loc(kernel_3845637904882483465))    
    call dict_add(dict_kernels, 3845637904882483465_8, kernel_3845637904882483465)
    !call dict_add(dict_kernels, 4814849100277162046_8, loc(kernel_4814849100277162046))    
    call dict_add(dict_kernels, 4814849100277162046_8, kernel_4814849100277162046)
    !call dict_add(dict_kernels, 4303050854118975013_8, loc(kernel_4303050854118975013))    
    call dict_add(dict_kernels, 4303050854118975013_8, kernel_4303050854118975013)
    !call dict_add(dict_kernels, 5391642689485582160_8, loc(kernel_5391642689485582160))    
    call dict_add(dict_kernels, 5391642689485582160_8, kernel_5391642689485582160)
    !call dict_add(dict_kernels, 8831349452088620100_8, loc(kernel_8831349452088620100))    
    call dict_add(dict_kernels, 8831349452088620100_8, kernel_8831349452088620100)
    !call dict_add(dict_kernels, 5222561026473440589_8, loc(kernel_5222561026473440589))    
    call dict_add(dict_kernels, 5222561026473440589_8, kernel_5222561026473440589)
    
    !call disp(dict_kernels)    
    ierr = 0
  end subroutine


  !> function X1*(-((A1)*(A2)))
  subroutine kernel_5319368613274967336(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)


    ! call tic()
    ! res(:,:,:) = ops_beta(1)*(-((x1(i,j,k))*(x2(i,j,k))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ops_beta(1)*(-((x1(i,j,k))*(x2(i,j,k))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    
  end subroutine

  !> function ((((X1*A1+Y1)-(X2*A2))+(X3*A3))-(X4*A4))+((X5*A5)*(A6))
  subroutine kernel_4639977022766068446(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)


    ! call tic()
    ! res(:,:,:) = ((((ops_beta(1)*x1(i,j,k)+ops_alpha(1))-(ops_beta(2)*x2(i,j,k)))+(ops_beta(3)*x3(i,j,k)))-(ops_beta(4)*x4(i,j,k)))+((ops_beta(5)*x5(i,j,k))*(x6(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((((ops_beta(1)*x1(i,j,k)+ops_alpha(1))-(ops_beta(2)*x2(i,j,k)))+(ops_beta(3)*x3(i,j,k)))-(ops_beta(4)*&
                 &x4(i,j,k)))+((ops_beta(5)*x5(i,j,k))*(x6(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    
  end subroutine

  !> function (((A1)+(((((X1*A2+Y1)+(X2*A3))-(X3*A4))+(X4*A5))*(A6)))+(((X5*A7+Y2)-(X6*A8))*((abs(A9))**N1)))+((X7*A10)*(A11))
  subroutine kernel_5104477621180273533(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9
    real(8), pointer :: x10(:,:,:)
    type(array), pointer :: B10
    real(8), pointer :: x11(:,:,:)
    type(array), pointer :: B11

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)
    B10 => null()    
    if(is_equal(A%dmda, B(10)%ptr%dmda)) then
       x10 => B(10)%ptr%data       
    else
       print*, "A is not equal to B(10)"
       call transfer(B10, A, B(10)%ptr)
       x10 => B10%data
    end if
    !print*, "x10 = ", lbound(x10), ubound(x10)
    B11 => null()    
    if(is_equal(A%dmda, B(11)%ptr%dmda)) then
       x11 => B(11)%ptr%data       
    else
       print*, "A is not equal to B(11)"
       call transfer(B11, A, B(11)%ptr)
       x11 => B11%data
    end if
    !print*, "x11 = ", lbound(x11), ubound(x11)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))+(((((ops_beta(1)*x2(i,j,k)+ops_alpha(1))+(ops_beta(2)*x3(i,j,k)))-(ops_beta(3)*x4(i,j,k)))+(ops_beta(4)*x5(i,j,k)))*(x6(i,j,k))))+(((ops_beta(5)*x7(i,j,k)+ops_alpha(2))-(ops_beta(6)*x8(i,j,k)))*((abs(x9(i,j,k)))**ops_args(1))))+((ops_beta(7)*x10(i,j,k))*(x11(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))+(((((ops_beta(1)*x2(i,j,k)+ops_alpha(1))+(ops_beta(2)*x3(i,j,k)))-(ops_beta(3)*x4(i,j,k))&
                 &)+(ops_beta(4)*x5(i,j,k)))*(x6(i,j,k))))+(((ops_beta(5)*x7(i,j,k)+ops_alpha(2))-(ops_beta(6)*x8(i,j,k)))*((abs(x9&
                 &(i,j,k)))**ops_args(1))))+((ops_beta(7)*x10(i,j,k))*(x11(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    call release_ptr(B10)
    call release_ptr(B11)
    
  end subroutine

  !> function (((X1*A1+Y1)+(X2*A2))-(X3*A3))+(X4*A4+Y2)
  subroutine kernel_3052836676995476825(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = (((ops_beta(1)*x1(i,j,k)+ops_alpha(1))+(ops_beta(2)*x2(i,j,k)))-(ops_beta(3)*x3(i,j,k)))+(ops_beta(4)*x4(i,j,k)+ops_alpha(2))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((ops_beta(1)*x1(i,j,k)+ops_alpha(1))+(ops_beta(2)*x2(i,j,k)))-(ops_beta(3)*x3(i,j,k)))+(ops_beta(4)*x&
                 &4(i,j,k)+ops_alpha(2))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function (A1)/((A2)*(A3))
  subroutine kernel_3254872525222050683(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))/((x2(i,j,k))*(x3(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))/((x2(i,j,k))*(x3(i,j,k)))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function (A1)+((X1*A2)*(X2*A3+Y1))
  subroutine kernel_4066808250212656659(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))+((ops_beta(1)*x2(i,j,k))*(ops_beta(2)*x3(i,j,k)+ops_alpha(1)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))+((ops_beta(1)*x2(i,j,k))*(ops_beta(2)*x3(i,j,k)+ops_alpha(1)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function (A1)*(sqrt(A2))
  subroutine kernel_3262871250189963503(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))*(sqrt(x2(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))*(sqrt(x2(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    
  end subroutine

  !> function (X2*(1.0/log(X1*A1)))**N1
  subroutine kernel_3384383430676826436(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)


    ! call tic()
    ! res(:,:,:) = (ops_beta(2)*(1.0/log(ops_beta(1)*x1(i,j,k))))**ops_args(1)
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (ops_beta(2)*(1.0/log(ops_beta(1)*x1(i,j,k))))**ops_args(1)
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    
  end subroutine

  !> function (X1*A1)+((A2)*(A3))
  subroutine kernel_1919310359753992010(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (ops_beta(1)*x1(i,j,k))+((x2(i,j,k))*(x3(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (ops_beta(1)*x1(i,j,k))+((x2(i,j,k))*(x3(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((A1)*(A2))+(X1*A3)
  subroutine kernel_5522439168916705770(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))*(x2(i,j,k)))+(ops_beta(1)*x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))*(x2(i,j,k)))+(ops_beta(1)*x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function -((A1)*(A2+Y1))
  subroutine kernel_3692924362984551671(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)


    ! call tic()
    ! res(:,:,:) = -((x1(i,j,k))*(x2(i,j,k)+ops_alpha(1)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = -((x1(i,j,k))*(x2(i,j,k)+ops_alpha(1)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    
  end subroutine

  !> function ((A1)*(A2))*((A3)+(A4))
  subroutine kernel_7913083383874115848(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))*(x2(i,j,k)))*((x3(i,j,k))+(x4(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))*(x2(i,j,k)))*((x3(i,j,k))+(x4(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function ((A1)+(A2))-(A3)
  subroutine kernel_646399903993819770(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))+(x2(i,j,k)))-(x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))+(x2(i,j,k)))-(x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((A1)*(A2))-((A3)*(A4))
  subroutine kernel_7761457064372070790(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))*(x2(i,j,k)))-((x3(i,j,k))*(x4(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))*(x2(i,j,k)))-((x3(i,j,k))*(x4(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function ((A1)*(A2))-(A3)
  subroutine kernel_644868325008555321(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))*(x2(i,j,k)))-(x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))*(x2(i,j,k)))-(x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((A1)+(A2))+(A3)
  subroutine kernel_646399903991447928(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))+(x2(i,j,k)))+(x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))+(x2(i,j,k)))+(x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((X1*A1)*(A2))*(A3)
  subroutine kernel_6766139683933094313(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((ops_beta(1)*x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((ops_beta(1)*x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function (-(((A1)*(A2))*(A3)))+(((A4)*(A5))*(A6))
  subroutine kernel_108513903478386821(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)


    ! call tic()
    ! res(:,:,:) = (-(((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k))))+(((x4(i,j,k))*(x5(i,j,k)))*(x6(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (-(((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k))))+(((x4(i,j,k))*(x5(i,j,k)))*(x6(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    
  end subroutine

  !> function ((X1*A1)*(A2))*(sqrt((((A3)**N1)+((A4)**N2))+(X2*(((A5)+(A6))**N3))))
  subroutine kernel_2381965364436835615(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)


    ! call tic()
    ! res(:,:,:) = ((ops_beta(1)*x1(i,j,k))*(x2(i,j,k)))*(sqrt((((x3(i,j,k))**ops_args(1))+((x4(i,j,k))**ops_args(2)))+(ops_beta(2)*(((x5(i,j,k))+(x6(i,j,k)))**ops_args(3)))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((ops_beta(1)*x1(i,j,k))*(x2(i,j,k)))*(sqrt((((x3(i,j,k))**ops_args(1))+((x4(i,j,k))**ops_args(2)))+(op&
                 &s_beta(2)*(((x5(i,j,k))+(x6(i,j,k)))**ops_args(3)))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    
  end subroutine

  !> function (A1)-(X1*(((A2)+(A3))-(A4)))
  subroutine kernel_5084845002459917026(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))-(ops_beta(1)*(((x2(i,j,k))+(x3(i,j,k)))-(x4(i,j,k))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))-(ops_beta(1)*(((x2(i,j,k))+(x3(i,j,k)))-(x4(i,j,k))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function (((A1)*(A2))-(X4*((((((A3)+(A4))-(A5))+((X1*A6)*(((X2*A7)+(X3*((A8)+(A9))))+(A10))))+(A11))+((A12)-(A13)))))/(A14)
  subroutine kernel_5404229555782944892(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9
    real(8), pointer :: x10(:,:,:)
    type(array), pointer :: B10
    real(8), pointer :: x11(:,:,:)
    type(array), pointer :: B11
    real(8), pointer :: x12(:,:,:)
    type(array), pointer :: B12
    real(8), pointer :: x13(:,:,:)
    type(array), pointer :: B13
    real(8), pointer :: x14(:,:,:)
    type(array), pointer :: B14

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)
    B10 => null()    
    if(is_equal(A%dmda, B(10)%ptr%dmda)) then
       x10 => B(10)%ptr%data       
    else
       print*, "A is not equal to B(10)"
       call transfer(B10, A, B(10)%ptr)
       x10 => B10%data
    end if
    !print*, "x10 = ", lbound(x10), ubound(x10)
    B11 => null()    
    if(is_equal(A%dmda, B(11)%ptr%dmda)) then
       x11 => B(11)%ptr%data       
    else
       print*, "A is not equal to B(11)"
       call transfer(B11, A, B(11)%ptr)
       x11 => B11%data
    end if
    !print*, "x11 = ", lbound(x11), ubound(x11)
    B12 => null()    
    if(is_equal(A%dmda, B(12)%ptr%dmda)) then
       x12 => B(12)%ptr%data       
    else
       print*, "A is not equal to B(12)"
       call transfer(B12, A, B(12)%ptr)
       x12 => B12%data
    end if
    !print*, "x12 = ", lbound(x12), ubound(x12)
    B13 => null()    
    if(is_equal(A%dmda, B(13)%ptr%dmda)) then
       x13 => B(13)%ptr%data       
    else
       print*, "A is not equal to B(13)"
       call transfer(B13, A, B(13)%ptr)
       x13 => B13%data
    end if
    !print*, "x13 = ", lbound(x13), ubound(x13)
    B14 => null()    
    if(is_equal(A%dmda, B(14)%ptr%dmda)) then
       x14 => B(14)%ptr%data       
    else
       print*, "A is not equal to B(14)"
       call transfer(B14, A, B(14)%ptr)
       x14 => B14%data
    end if
    !print*, "x14 = ", lbound(x14), ubound(x14)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*((((((x3(i,j,k))+(x4(i,j,k)))-(x5(i,j,k)))+((ops_beta(1)*x6(i,j,k))*(((ops_beta(2)*x7(i,j,k))+(ops_beta(3)*((x8(i,j,k))+(x9(i,j,k)))))+(x10(i,j,k)))))+(x11(i,j,k)))+((x12(i,j,k))-(x13(i,j,k))))))/(x14(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*((((((x3(i,j,k))+(x4(i,j,k)))-(x5(i,j,k)))+((ops_beta(1)*x6(i,j&
                 &,k))*(((ops_beta(2)*x7(i,j,k))+(ops_beta(3)*((x8(i,j,k))+(x9(i,j,k)))))+(x10(i,j,k)))))+(x11(i,j,k)))+((x12(i,j,k&
                 &))-(x13(i,j,k))))))/(x14(i,j,k))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    call release_ptr(B10)
    call release_ptr(B11)
    call release_ptr(B12)
    call release_ptr(B13)
    call release_ptr(B14)
    
  end subroutine

  !> function ((A1)*(A2))*(A3)
  subroutine kernel_644868325004997558(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function X3*((A1)+((X2*(sqrt(X1*(1.0/A2))))*((A3)-(A4))))
  subroutine kernel_7381972720285064861(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = ops_beta(3)*((x1(i,j,k))+((ops_beta(2)*(sqrt(ops_beta(1)*(1.0/x2(i,j,k)))))*((x3(i,j,k))-(x4(i,j,k)))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ops_beta(3)*((x1(i,j,k))+((ops_beta(2)*(sqrt(ops_beta(1)*(1.0/x2(i,j,k)))))*((x3(i,j,k))-(x4(i,j,k)))))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function X3*((A1)-((X2*(sqrt(X1*(1.0/A2))))*((A3)-(A4))))
  subroutine kernel_8205386001274835749(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = ops_beta(3)*((x1(i,j,k))-((ops_beta(2)*(sqrt(ops_beta(1)*(1.0/x2(i,j,k)))))*((x3(i,j,k))-(x4(i,j,k)))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ops_beta(3)*((x1(i,j,k))-((ops_beta(2)*(sqrt(ops_beta(1)*(1.0/x2(i,j,k)))))*((x3(i,j,k))-(x4(i,j,k)))))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function (((A1)*(A2))-(X4*((((((A3)+(A4))+(A5))+((X1*A6)*(((X2*A7)+(X3*((A8)+(A9))))+(A10))))+(A11))+((A12)-(A13)))))/(A14)
  subroutine kernel_5507350668782427326(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9
    real(8), pointer :: x10(:,:,:)
    type(array), pointer :: B10
    real(8), pointer :: x11(:,:,:)
    type(array), pointer :: B11
    real(8), pointer :: x12(:,:,:)
    type(array), pointer :: B12
    real(8), pointer :: x13(:,:,:)
    type(array), pointer :: B13
    real(8), pointer :: x14(:,:,:)
    type(array), pointer :: B14

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)
    B10 => null()    
    if(is_equal(A%dmda, B(10)%ptr%dmda)) then
       x10 => B(10)%ptr%data       
    else
       print*, "A is not equal to B(10)"
       call transfer(B10, A, B(10)%ptr)
       x10 => B10%data
    end if
    !print*, "x10 = ", lbound(x10), ubound(x10)
    B11 => null()    
    if(is_equal(A%dmda, B(11)%ptr%dmda)) then
       x11 => B(11)%ptr%data       
    else
       print*, "A is not equal to B(11)"
       call transfer(B11, A, B(11)%ptr)
       x11 => B11%data
    end if
    !print*, "x11 = ", lbound(x11), ubound(x11)
    B12 => null()    
    if(is_equal(A%dmda, B(12)%ptr%dmda)) then
       x12 => B(12)%ptr%data       
    else
       print*, "A is not equal to B(12)"
       call transfer(B12, A, B(12)%ptr)
       x12 => B12%data
    end if
    !print*, "x12 = ", lbound(x12), ubound(x12)
    B13 => null()    
    if(is_equal(A%dmda, B(13)%ptr%dmda)) then
       x13 => B(13)%ptr%data       
    else
       print*, "A is not equal to B(13)"
       call transfer(B13, A, B(13)%ptr)
       x13 => B13%data
    end if
    !print*, "x13 = ", lbound(x13), ubound(x13)
    B14 => null()    
    if(is_equal(A%dmda, B(14)%ptr%dmda)) then
       x14 => B(14)%ptr%data       
    else
       print*, "A is not equal to B(14)"
       call transfer(B14, A, B(14)%ptr)
       x14 => B14%data
    end if
    !print*, "x14 = ", lbound(x14), ubound(x14)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*((((((x3(i,j,k))+(x4(i,j,k)))+(x5(i,j,k)))+((ops_beta(1)*x6(i,j,k))*(((ops_beta(2)*x7(i,j,k))+(ops_beta(3)*((x8(i,j,k))+(x9(i,j,k)))))+(x10(i,j,k)))))+(x11(i,j,k)))+((x12(i,j,k))-(x13(i,j,k))))))/(x14(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*((((((x3(i,j,k))+(x4(i,j,k)))+(x5(i,j,k)))+((ops_beta(1)*x6(i,j&
                 &,k))*(((ops_beta(2)*x7(i,j,k))+(ops_beta(3)*((x8(i,j,k))+(x9(i,j,k)))))+(x10(i,j,k)))))+(x11(i,j,k)))+((x12(i,j,k&
                 &))-(x13(i,j,k))))))/(x14(i,j,k))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    call release_ptr(B10)
    call release_ptr(B11)
    call release_ptr(B12)
    call release_ptr(B13)
    call release_ptr(B14)
    
  end subroutine

  !> function (A1)+(X2*(((A2)-(X1*A3))+(A4)))
  subroutine kernel_3303401064826175212(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))+(ops_beta(2)*(((x2(i,j,k))-(ops_beta(1)*x3(i,j,k)))+(x4(i,j,k))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))+(ops_beta(2)*(((x2(i,j,k))-(ops_beta(1)*x3(i,j,k)))+(x4(i,j,k))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function (A1)+(X2*((X1*A2)*(A3)))
  subroutine kernel_2281688410210494257(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))+(ops_beta(2)*((ops_beta(1)*x2(i,j,k))*(x3(i,j,k))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))+(ops_beta(2)*((ops_beta(1)*x2(i,j,k))*(x3(i,j,k))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((A1)*(A2))-(((X1*A3)*(A4))*(A5))
  subroutine kernel_259121922155312303(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))*(x2(i,j,k)))-(((ops_beta(1)*x3(i,j,k))*(x4(i,j,k)))*(x5(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))*(x2(i,j,k)))-(((ops_beta(1)*x3(i,j,k))*(x4(i,j,k)))*(x5(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    
  end subroutine

  !> function ((A1)*(A2))-(((A3)*(A4))*((A5)+(A6)))
  subroutine kernel_7595890524641150240(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))*(x2(i,j,k)))-(((x3(i,j,k))*(x4(i,j,k)))*((x5(i,j,k))+(x6(i,j,k))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))*(x2(i,j,k)))-(((x3(i,j,k))*(x4(i,j,k)))*((x5(i,j,k))+(x6(i,j,k))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    
  end subroutine

  !> function ((A1)+(X1*A2))*(A3)
  subroutine kernel_6383578275131078730(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))+(ops_beta(1)*x2(i,j,k)))*(x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))+(ops_beta(1)*x2(i,j,k)))*(x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((((A1)*(A2))*(A3))*(A4))*(A5)
  subroutine kernel_8877419281824194841(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)


    ! call tic()
    ! res(:,:,:) = ((((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k)))*(x4(i,j,k)))*(x5(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k)))*(x4(i,j,k)))*(x5(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    
  end subroutine

  !> function (X1*((A1)+(A2)))/(A3)
  subroutine kernel_1039351428668170720(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (ops_beta(1)*((x1(i,j,k))+(x2(i,j,k))))/(x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (ops_beta(1)*((x1(i,j,k))+(x2(i,j,k))))/(x3(i,j,k))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function (A1)*(((A2)+(A3))+(X1*((A4)-(A5))))
  subroutine kernel_1279377163763675870(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))*(((x2(i,j,k))+(x3(i,j,k)))+(ops_beta(1)*((x4(i,j,k))-(x5(i,j,k)))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))*(((x2(i,j,k))+(x3(i,j,k)))+(ops_beta(1)*((x4(i,j,k))-(x5(i,j,k)))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    
  end subroutine

  !> function (((A1)*(A2))-(X1*(((-(A3))+(A4))+(A5))))/(A6)
  subroutine kernel_8707219604264828824(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(1)*(((-(x3(i,j,k)))+(x4(i,j,k)))+(x5(i,j,k)))))/(x6(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(1)*(((-(x3(i,j,k)))+(x4(i,j,k)))+(x5(i,j,k)))))/(x6(i,j,k))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    
  end subroutine

  !> function sqrt(((A1)**N1)+((A2)**N2))
  subroutine kernel_4905021281077427515(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)


    ! call tic()
    ! res(:,:,:) = sqrt(((x1(i,j,k))**ops_args(1))+((x2(i,j,k))**ops_args(2)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = sqrt(((x1(i,j,k))**ops_args(1))+((x2(i,j,k))**ops_args(2)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    
  end subroutine

  !> function (((X1*A1+Y1)+(X2*A2+Y2))-(X3*((A3+Y3)**N1)))+(X4*A4+Y4)
  subroutine kernel_1819205572515252931(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = (((ops_beta(1)*x1(i,j,k)+ops_alpha(1))+(ops_beta(2)*x2(i,j,k)+ops_alpha(2)))-(ops_beta(3)*((x3(i,j,k)+ops_alpha(3))**ops_args(1))))+(ops_beta(4)*x4(i,j,k)+ops_alpha(4))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((ops_beta(1)*x1(i,j,k)+ops_alpha(1))+(ops_beta(2)*x2(i,j,k)+ops_alpha(2)))-(ops_beta(3)*((x3(i,j,k)+o&
                 &ps_alpha(3))**ops_args(1))))+(ops_beta(4)*x4(i,j,k)+ops_alpha(4))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function (A1)/(sqrt((X2*((X1*A2)/(A3))+Y1)*(X4*((X3*A4)/((A5)**N1))+Y2)))
  subroutine kernel_5940693024572846470(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))/(sqrt((ops_beta(2)*((ops_beta(1)*x2(i,j,k))/(x3(i,j,k)))+ops_alpha(1))*(ops_beta(4)*((ops_beta(3)*x4(i,j,k))/((x5(i,j,k))**ops_args(1)))+ops_alpha(2))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))/(sqrt((ops_beta(2)*((ops_beta(1)*x2(i,j,k))/(x3(i,j,k)))+ops_alpha(1))*(ops_beta(4)*((ops_b&
                 &eta(3)*x4(i,j,k))/((x5(i,j,k))**ops_args(1)))+ops_alpha(2))))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    
  end subroutine

  !> function (-((X1*A1)/(A2)))+(X2*(1.0/A3))
  subroutine kernel_3864818436564201008(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (-((ops_beta(1)*x1(i,j,k))/(x2(i,j,k))))+(ops_beta(2)*(1.0/x3(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (-((ops_beta(1)*x1(i,j,k))/(x2(i,j,k))))+(ops_beta(2)*(1.0/x3(i,j,k)))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((((X1*A1)*(((A2)**N1)+((A3)**N2)))/((A4)**N3))-((X2*A5)*(A6)))+((A7)*(A8))
  subroutine kernel_1820361344622479259(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)


    ! call tic()
    ! res(:,:,:) = ((((ops_beta(1)*x1(i,j,k))*(((x2(i,j,k))**ops_args(1))+((x3(i,j,k))**ops_args(2))))/((x4(i,j,k))**ops_args(3)))-((ops_beta(2)*x5(i,j,k))*(x6(i,j,k))))+((x7(i,j,k))*(x8(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((((ops_beta(1)*x1(i,j,k))*(((x2(i,j,k))**ops_args(1))+((x3(i,j,k))**ops_args(2))))/((x4(i,j,k))**ops_a&
                 &rgs(3)))-((ops_beta(2)*x5(i,j,k))*(x6(i,j,k))))+((x7(i,j,k))*(x8(i,j,k)))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    
  end subroutine

  !> function (sqrt(A1))/(X1*A2+Y1)
  subroutine kernel_2483251290183065470(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)


    ! call tic()
    ! res(:,:,:) = (sqrt(x1(i,j,k)))/(ops_beta(1)*x2(i,j,k)+ops_alpha(1))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (sqrt(x1(i,j,k)))/(ops_beta(1)*x2(i,j,k)+ops_alpha(1))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    
  end subroutine

  !> function (X1*A1+Y1)/(X3*((X2*(1.0/A2)+Y2)*(A3))+Y3)
  subroutine kernel_8238309252147448007(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (ops_beta(1)*x1(i,j,k)+ops_alpha(1))/(ops_beta(3)*((ops_beta(2)*(1.0/x2(i,j,k))+ops_alpha(2))*(x3(i,j,k)))+ops_alpha(3))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (ops_beta(1)*x1(i,j,k)+ops_alpha(1))/(ops_beta(3)*((ops_beta(2)*(1.0/x2(i,j,k))+ops_alpha(2))*(x3(i,j,k&
                 &)))+ops_alpha(3))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((X1*A1+Y1)+((X2*A2)*(A3)))/(X3*A4+Y2)
  subroutine kernel_445431549541144228(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = ((ops_beta(1)*x1(i,j,k)+ops_alpha(1))+((ops_beta(2)*x2(i,j,k))*(x3(i,j,k))))/(ops_beta(3)*x4(i,j,k)+ops_alpha(2))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((ops_beta(1)*x1(i,j,k)+ops_alpha(1))+((ops_beta(2)*x2(i,j,k))*(x3(i,j,k))))/(ops_beta(3)*x4(i,j,k)+ops&
                 &_alpha(2))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function (((A1)*(A2))-(X1*(((A3)+(A4))-(A5))))/(A6)
  subroutine kernel_6348923055873632844(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(1)*(((x3(i,j,k))+(x4(i,j,k)))-(x5(i,j,k)))))/(x6(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(1)*(((x3(i,j,k))+(x4(i,j,k)))-(x5(i,j,k)))))/(x6(i,j,k))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    
  end subroutine

  !> function (((A1)*(A2))*(A3))-(((X1*((A4)*(A5)))*(A6))*(A7))
  subroutine kernel_4164740083217831399(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k)))-(((ops_beta(1)*((x4(i,j,k))*(x5(i,j,k))))*(x6(i,j,k)))*(x7(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))*(x2(i,j,k)))*(x3(i,j,k)))-(((ops_beta(1)*((x4(i,j,k))*(x5(i,j,k))))*(x6(i,j,k)))*(x7(i,j,&
                 &k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    
  end subroutine

  !> function (A1)-(X3*((((X1*((A2)-(abs(A3))))*((A4)-(A5)))/(A6))+((X2*((A7)+(abs(A8))))*(A9))))
  subroutine kernel_7821363202111267096(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))-(ops_beta(3)*((((ops_beta(1)*((x2(i,j,k))-(abs(x3(i,j,k)))))*((x4(i,j,k))-(x5(i,j,k))))/(x6(i,j,k)))+((ops_beta(2)*((x7(i,j,k))+(abs(x8(i,j,k)))))*(x9(i,j,k)))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))-(ops_beta(3)*((((ops_beta(1)*((x2(i,j,k))-(abs(x3(i,j,k)))))*((x4(i,j,k))-(x5(i,j,k))))/(x6&
                 &(i,j,k)))+((ops_beta(2)*((x7(i,j,k))+(abs(x8(i,j,k)))))*(x9(i,j,k)))))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    
  end subroutine

  !> function (A1)+(((((((X1*((A2)+(abs(A3))))/(A4))*(A5))/(A6))*(A7))*(A8))/(A9))
  subroutine kernel_7381297727881908921(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))+(((((((ops_beta(1)*((x2(i,j,k))+(abs(x3(i,j,k)))))/(x4(i,j,k)))*(x5(i,j,k)))/(x6(i,j,k)))*(x7(i,j,k)))*(x8(i,j,k)))/(x9(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))+(((((((ops_beta(1)*((x2(i,j,k))+(abs(x3(i,j,k)))))/(x4(i,j,k)))*(x5(i,j,k)))/(x6(i,j,k)))*(&
                 &x7(i,j,k)))*(x8(i,j,k)))/(x9(i,j,k)))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    
  end subroutine

  !> function (A1)-(X3*((((X1*((A2)+(abs(A3))))*((A4)-(A5)))/(A6))+((((X2*((A7)-(abs(A8))))*(A9))*(A10))/(A11))))
  subroutine kernel_563110658223512454(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9
    real(8), pointer :: x10(:,:,:)
    type(array), pointer :: B10
    real(8), pointer :: x11(:,:,:)
    type(array), pointer :: B11

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)
    B10 => null()    
    if(is_equal(A%dmda, B(10)%ptr%dmda)) then
       x10 => B(10)%ptr%data       
    else
       print*, "A is not equal to B(10)"
       call transfer(B10, A, B(10)%ptr)
       x10 => B10%data
    end if
    !print*, "x10 = ", lbound(x10), ubound(x10)
    B11 => null()    
    if(is_equal(A%dmda, B(11)%ptr%dmda)) then
       x11 => B(11)%ptr%data       
    else
       print*, "A is not equal to B(11)"
       call transfer(B11, A, B(11)%ptr)
       x11 => B11%data
    end if
    !print*, "x11 = ", lbound(x11), ubound(x11)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))-(ops_beta(3)*((((ops_beta(1)*((x2(i,j,k))+(abs(x3(i,j,k)))))*((x4(i,j,k))-(x5(i,j,k))))/(x6(i,j,k)))+((((ops_beta(2)*((x7(i,j,k))-(abs(x8(i,j,k)))))*(x9(i,j,k)))*(x10(i,j,k)))/(x11(i,j,k)))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))-(ops_beta(3)*((((ops_beta(1)*((x2(i,j,k))+(abs(x3(i,j,k)))))*((x4(i,j,k))-(x5(i,j,k))))/(x6&
                 &(i,j,k)))+((((ops_beta(2)*((x7(i,j,k))-(abs(x8(i,j,k)))))*(x9(i,j,k)))*(x10(i,j,k)))/(x11(i,j,k)))))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    call release_ptr(B10)
    call release_ptr(B11)
    
  end subroutine

  !> function (A1)+(((((((X1*((A2)-(abs(A3))))/(A4))*(A5))/(A6))*(A7))*(A8))/(A9))
  subroutine kernel_3845637904882483465(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))+(((((((ops_beta(1)*((x2(i,j,k))-(abs(x3(i,j,k)))))/(x4(i,j,k)))*(x5(i,j,k)))/(x6(i,j,k)))*(x7(i,j,k)))*(x8(i,j,k)))/(x9(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))+(((((((ops_beta(1)*((x2(i,j,k))-(abs(x3(i,j,k)))))/(x4(i,j,k)))*(x5(i,j,k)))/(x6(i,j,k)))*(&
                 &x7(i,j,k)))*(x8(i,j,k)))/(x9(i,j,k)))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    
  end subroutine

  !> function (((A1)*(A2))-(X4*(((((A3)+(A4))-(A5))+(X3*((X1*A6)*((A7)+(X2*A8)))))-(A9))))/(A10)
  subroutine kernel_4814849100277162046(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9
    real(8), pointer :: x10(:,:,:)
    type(array), pointer :: B10

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)
    B10 => null()    
    if(is_equal(A%dmda, B(10)%ptr%dmda)) then
       x10 => B(10)%ptr%data       
    else
       print*, "A is not equal to B(10)"
       call transfer(B10, A, B(10)%ptr)
       x10 => B10%data
    end if
    !print*, "x10 = ", lbound(x10), ubound(x10)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*(((((x3(i,j,k))+(x4(i,j,k)))-(x5(i,j,k)))+(ops_beta(3)*((ops_beta(1)*x6(i,j,k))*((x7(i,j,k))+(ops_beta(2)*x8(i,j,k))))))-(x9(i,j,k)))))/(x10(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*(((((x3(i,j,k))+(x4(i,j,k)))-(x5(i,j,k)))+(ops_beta(3)*((ops_be&
                 &ta(1)*x6(i,j,k))*((x7(i,j,k))+(ops_beta(2)*x8(i,j,k))))))-(x9(i,j,k)))))/(x10(i,j,k))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    call release_ptr(B10)
    
  end subroutine

  !> function (A1)*(sqrt(((A2)**N1)+((A3)**N2)))
  subroutine kernel_4303050854118975013(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = (x1(i,j,k))*(sqrt(((x2(i,j,k))**ops_args(1))+((x3(i,j,k))**ops_args(2))))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (x1(i,j,k))*(sqrt(((x2(i,j,k))**ops_args(1))+((x3(i,j,k))**ops_args(2))))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine

  !> function ((sqrt(X1*A1))*(A2))+((X3*(sqrt(X2*A3))+Y1)*(A4))
  subroutine kernel_5391642689485582160(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)


    ! call tic()
    ! res(:,:,:) = ((sqrt(ops_beta(1)*x1(i,j,k)))*(x2(i,j,k)))+((ops_beta(3)*(sqrt(ops_beta(2)*x3(i,j,k)))+ops_alpha(1))*(x4(i,j,k)))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((sqrt(ops_beta(1)*x1(i,j,k)))*(x2(i,j,k)))+((ops_beta(3)*(sqrt(ops_beta(2)*x3(i,j,k)))+ops_alpha(1))*(&
                 &x4(i,j,k)))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    
  end subroutine

  !> function (((A1)*(A2))-(X4*(((((A3)+(A4))+(A5))+(X3*((X1*A6)*((A7)+(X2*A8)))))-(A9))))/(A10)
  subroutine kernel_8831349452088620100(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3
    real(8), pointer :: x4(:,:,:)
    type(array), pointer :: B4
    real(8), pointer :: x5(:,:,:)
    type(array), pointer :: B5
    real(8), pointer :: x6(:,:,:)
    type(array), pointer :: B6
    real(8), pointer :: x7(:,:,:)
    type(array), pointer :: B7
    real(8), pointer :: x8(:,:,:)
    type(array), pointer :: B8
    real(8), pointer :: x9(:,:,:)
    type(array), pointer :: B9
    real(8), pointer :: x10(:,:,:)
    type(array), pointer :: B10

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)
    B4 => null()    
    if(is_equal(A%dmda, B(4)%ptr%dmda)) then
       x4 => B(4)%ptr%data       
    else
       print*, "A is not equal to B(4)"
       call transfer(B4, A, B(4)%ptr)
       x4 => B4%data
    end if
    !print*, "x4 = ", lbound(x4), ubound(x4)
    B5 => null()    
    if(is_equal(A%dmda, B(5)%ptr%dmda)) then
       x5 => B(5)%ptr%data       
    else
       print*, "A is not equal to B(5)"
       call transfer(B5, A, B(5)%ptr)
       x5 => B5%data
    end if
    !print*, "x5 = ", lbound(x5), ubound(x5)
    B6 => null()    
    if(is_equal(A%dmda, B(6)%ptr%dmda)) then
       x6 => B(6)%ptr%data       
    else
       print*, "A is not equal to B(6)"
       call transfer(B6, A, B(6)%ptr)
       x6 => B6%data
    end if
    !print*, "x6 = ", lbound(x6), ubound(x6)
    B7 => null()    
    if(is_equal(A%dmda, B(7)%ptr%dmda)) then
       x7 => B(7)%ptr%data       
    else
       print*, "A is not equal to B(7)"
       call transfer(B7, A, B(7)%ptr)
       x7 => B7%data
    end if
    !print*, "x7 = ", lbound(x7), ubound(x7)
    B8 => null()    
    if(is_equal(A%dmda, B(8)%ptr%dmda)) then
       x8 => B(8)%ptr%data       
    else
       print*, "A is not equal to B(8)"
       call transfer(B8, A, B(8)%ptr)
       x8 => B8%data
    end if
    !print*, "x8 = ", lbound(x8), ubound(x8)
    B9 => null()    
    if(is_equal(A%dmda, B(9)%ptr%dmda)) then
       x9 => B(9)%ptr%data       
    else
       print*, "A is not equal to B(9)"
       call transfer(B9, A, B(9)%ptr)
       x9 => B9%data
    end if
    !print*, "x9 = ", lbound(x9), ubound(x9)
    B10 => null()    
    if(is_equal(A%dmda, B(10)%ptr%dmda)) then
       x10 => B(10)%ptr%data       
    else
       print*, "A is not equal to B(10)"
       call transfer(B10, A, B(10)%ptr)
       x10 => B10%data
    end if
    !print*, "x10 = ", lbound(x10), ubound(x10)


    ! call tic()
    ! res(:,:,:) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*(((((x3(i,j,k))+(x4(i,j,k)))+(x5(i,j,k)))+(ops_beta(3)*((ops_beta(1)*x6(i,j,k))*((x7(i,j,k))+(ops_beta(2)*x8(i,j,k))))))-(x9(i,j,k)))))/(x10(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = (((x1(i,j,k))*(x2(i,j,k)))-(ops_beta(4)*(((((x3(i,j,k))+(x4(i,j,k)))+(x5(i,j,k)))+(ops_beta(3)*((ops_be&
                 &ta(1)*x6(i,j,k))*((x7(i,j,k))+(ops_beta(2)*x8(i,j,k))))))-(x9(i,j,k)))))/(x10(i,j,k))
          end do
       end do
    end do

    do k = zs,ze
       do j = ys, ye
          do i = xs, xe
             if(i == xs .or. i == xe .or. &
                  j == ys .or. j == ye .or. &
                  k == zs .or. k == ze) then
                if(isnan(res(i, j, k)) .or. &
                     res(i,j,k)==inf) &
                  res(i,j,k) = 0
             endif
          end do
       end do
    end do
    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    call release_ptr(B4)
    call release_ptr(B5)
    call release_ptr(B6)
    call release_ptr(B7)
    call release_ptr(B8)
    call release_ptr(B9)
    call release_ptr(B10)
    
  end subroutine

  !> function ((A1)+(A2))-(X1*A3)
  subroutine kernel_5222561026473440589(A, B, ops_alpha, &
       ops_beta, ops_args, ops_num, args_num,  A_alpha, A_beta, ierr)
    implicit none
    type(array), intent(inout) :: A
    type(array_ptr), intent(in) :: B(:)
    real(8), intent(in) :: ops_alpha(:), A_alpha, A_beta
    real(8), intent(in) :: ops_beta(:)
    integer, intent(in)  :: args_num, ops_num    
    real(8), intent(in) :: ops_args(:)
    real(8), pointer :: res(:,:,:)
    integer :: xs, xe, ys, ye, zs, ze, i,j,k
    integer, intent(out) :: ierr
    real(8) :: t1, t2, inf
    real(8), pointer :: x1(:,:,:)
    type(array), pointer :: B1
    real(8), pointer :: x2(:,:,:)
    type(array), pointer :: B2
    real(8), pointer :: x3(:,:,:)
    type(array), pointer :: B3

    t1 = 1.0
    t2 = 0.0
    inf = t1 / t2
    
    res => A%data
    ! print*, "res = ", lbound(res), ubound(res)
    !print*, "ops_beta", ops_beta
    !print*, "ops_alpha", ops_alpha
    ! if(size(ops_args) > 0) then
    !    print*, "ops_args=", ops_args
    ! end if
    
    call set_corners(lbox(A), xs,ys,zs,xe,ye,ze)

    B1 => null()    
    if(is_equal(A%dmda, B(1)%ptr%dmda)) then
       x1 => B(1)%ptr%data       
    else
       print*, "A is not equal to B(1)"
       call transfer(B1, A, B(1)%ptr)
       x1 => B1%data
    end if
    !print*, "x1 = ", lbound(x1), ubound(x1)
    B2 => null()    
    if(is_equal(A%dmda, B(2)%ptr%dmda)) then
       x2 => B(2)%ptr%data       
    else
       print*, "A is not equal to B(2)"
       call transfer(B2, A, B(2)%ptr)
       x2 => B2%data
    end if
    !print*, "x2 = ", lbound(x2), ubound(x2)
    B3 => null()    
    if(is_equal(A%dmda, B(3)%ptr%dmda)) then
       x3 => B(3)%ptr%data       
    else
       print*, "A is not equal to B(3)"
       call transfer(B3, A, B(3)%ptr)
       x3 => B3%data
    end if
    !print*, "x3 = ", lbound(x3), ubound(x3)


    ! call tic()
    ! res(:,:,:) = ((x1(i,j,k))+(x2(i,j,k)))-(ops_beta(1)*x3(i,j,k))
    ! call toc()
    
    do k = zs,ze
       do j = ys, ye
          !dir$ simd 
          do i = xs, xe
             res(i, j, k) = ((x1(i,j,k))+(x2(i,j,k)))-(ops_beta(1)*x3(i,j,k))
          end do
       end do
    end do

    
    call release_ptr(B1)
    call release_ptr(B2)
    call release_ptr(B3)
    
  end subroutine
end module
