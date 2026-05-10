program fortran_ensight_test

   integer, parameter :: nnos = 8, nel=1
   real(4) :: x(nnos), y(nnos), z(nnos)
   character :: buffer*80
   
   x = 0.0
   y = 0.0
   z = 0.0
   
   x(2:3) = 1.0
   x(6:7) = 1.0
   
   y(3:4) = 1.0
   y(7:8) = 1.0
   
   z(5:8) = 1.0
     
     
   ! G E O M E T R Y
   !========================================================================

   open(9,file='foo.geo',form='unformatted')

   buffer = 'Fortran Binary'
   write (9) buffer
   buffer = 'description line 1'
   write (9) buffer
   buffer = 'description line 2'
   write (9) buffer
   buffer = 'node id given'
   write (9) buffer
   buffer = 'element id given'
   write (9) buffer
   buffer = 'part'
   write (9) buffer
   write (9) 1
   buffer = 'description line'
   write (9) buffer

   buffer = 'coordinates'
   write (9) buffer
   write (9) nnos
   write (9) (i, i=1,nnos)
   write (9) (x(i),i=1,nnos)
   write (9) (y(i),i=1,nnos)
   write (9) (z(i),i=1,nnos)
   
   buffer = 'hexa8'
   write (9) buffer
   write (9) nel
   write (9) (i, i=1,nel)
   write (9) (i, i=1,8)
   
   close(9)

   ! S C A L A R   F I E L D   P E R    N O D E 
   !========================================================================
   
   open(9,file='foo.scl',form='unformatted')

   buffer = 'scalar file'
   write(9) buffer
   buffer = 'part'
   write(9) buffer
   write(9) 1
   buffer = 'coordinates'
   write(9) buffer
   write(9) (sngl(i),i=1,nnos)
   close(9)
   
   close(9)

   ! V E C T O R   F I E L D   P E R    N O D E 
   !========================================================================

   open(9,file='foo.vec',form='unformatted')

   buffer = 'vector file'
   write(9) buffer
   buffer = 'part'
   write(9) buffer
   write(9) 1
   buffer = 'coordinates'
   write(9) buffer
   write(9) (x(i),i=1,nnos)
   write(9) (y(i),i=1,nnos)
   write(9) (z(i),i=1,nnos)

   close(9)

   ! C A S E    F I L E 
   !========================================================================
   
   open(9,file='foo.case')
   
   write(9,'(A)') 'FORMAT'
   write(9,'(A)') 'type: ensight gold'
   write(9,'(A)') 'GEOMETRY'
   write(9,'(A)') 'model: foo.geo'
   write(9,'(A)') 'VARIABLE'
   write(9,'(A)') 'scalar per node:   scalar   foo.scl'
   write(9,'(A)') 'vector per node:   vector   foo.vec'
   
   close(9)

end program
