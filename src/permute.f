      subroutine permute(class,pclass,numplt,tclass)
c
      integer numplt
      integer class(numplt)
      integer pclass(numplt)
c
c* local
c
      integer tclass(numplt)
      integer pool
c
      pool = numplt
      do 10 i=1,numplt
      tclass(i) = class(i)
   10 continue
c
      do 11 i=1,numplt
      index = rand()*(pool)+1
      pclass(i) = tclass(index)
      tclass(index) = tclass(pool)
      pool = pool - 1
   11 continue
c
      return
c
      end

