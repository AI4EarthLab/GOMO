#include "common.h"

subroutine bcond1()
  use openarray
  use variables
  use config
  implicit none
    
  call set(sub(elf,  1 , ':'), sub(elf, 2   , ':'))
  call set(sub(elf,  im, ':'), sub(elf, imm1, ':'))
  call set(sub(elf, ':',  1 ), sub(elf, ':' ,  2 ))
  call set(sub(elf, ':',  jm), sub(elf, ':' , jmm1))
  elf=elf*fsm

end subroutine bcond1
